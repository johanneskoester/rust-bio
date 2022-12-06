// Copyright 2021 Todd Morse.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Traits and utilities to read and write files in the FASTA and FASTQ format interchangably
//!
//! Files in the FASTA and FASTQ format have several common fields: ID, Sequence, and Description.
//! These utilities can be used to implement algorithms that only require these fields for files
//! that are in either format.
//!
//! This module serves two use cases:
//!
//! 1. Implementing functions that can be used generically with FASTA/FASTQ records. In this use
//! case the type may be known at compile time by the caller of your function.
//! 2. Processing data that may be either in the FASTA/FASTQ format. In this use case the type
//! cannot be known at compile time and you may or may not want to treat FASTA/FASTQ data
//! differently.
//!
//! # Generic Implementation Examples
//!
//! ## Common Statistics
//!
//! In this example, we implement a `count_bases` function that for supports both the FASTA and FASTQ
//! format.
//!
//! ```
//! use bio::io::{fasta, fastq, fastx};
//! use std::io;
//!
//! fn count_bases<T, E, I>(records: I) -> Result<usize, E>
//! where
//!     T: fastx::Record,
//!     E: std::error::Error,
//!     I: fastx::Records<T, E>,
//! {
//!     let mut nb_bases = 0;
//!     for result in records {
//!         let record = result?;
//!         nb_bases += record.seq().len();
//!     }
//!     Ok(nb_bases)
//! }
//!
//! let mut raw_reader = io::Cursor::new(
//!     b">id desc
//! ACTG
//! ",
//! );
//!
//! match fastx::get_kind_seek(&mut raw_reader) {
//!     Ok(fastx::Kind::FASTA) => {
//!         let mut reader = fasta::Reader::new(raw_reader);
//!         let nb_bases = count_bases(reader.records()).unwrap();
//!         println!("Number of bases: {}", nb_bases);
//!     }
//!     Ok(fastx::Kind::FASTQ) => {
//!         let mut reader = fastq::Reader::new(raw_reader);
//!         let nb_bases = count_bases(reader.records()).unwrap();
//!         println!("Number of bases: {}", nb_bases);
//!     }
//!     _ => println!("Encountered an error"),
//! }
//! ```
//!
//! ## Filtration
//!
//! In this example, we define an `at_least_n_bases` function that can filter FASTA or FASTQ
//! records based on their sequence lengths. It works seemlessly with `fasta::Record`s as if
//! it was implemented just for them. In a realistic scenario this function might be
//! defined in a library so callers could use it with both FASTA and FASTQ files as needed.
//!
//! ```
//! use bio::io::{fasta, fastq, fastx};
//! use std::io;
//! use std::io::BufReader;
//!
//! fn at_least_n_bases<T, E, I>(records: I, n: usize) -> impl Iterator<Item = Result<T, E>>
//! where
//!     T: fastx::Record,
//!     E: std::error::Error,
//!     I: fastx::Records<T, E>,
//! {
//!     records.filter(move |rr| match rr {
//!         Ok(r) => r.seq().len() > n,
//!         _ => true,
//!     })
//! }
//!
//! let mut reader = fasta::Reader::new(BufReader::new(io::stdin()));
//! let mut writer = fasta::Writer::new(io::stdout());
//!
//! for record in at_least_n_bases(reader.records(), 10) {
//!     writer.write_record(&record.unwrap());
//! }
//! ```
//! # Dynamic Type Examples
//!
//! If the type of a record is not known at compile time you can represent it with `fastx::EitherRecord`.
//! This type is an enum containing either a `fasta::Record` or a `fastq::Record`. There are also
//! utility functions defined on the enum so you can work with them without converting to the
//! underlying type.
//!
//! Note that using dynamic types like this incurs a slight performance cost around 5-10% when
//! compared to the typing based approach or parsing a single type directly.
//!
//! ## Parsing data of either the FASTA of FASTQ type
//!
//! ```
//! use bio::io::fastx::{EitherRecords, Record};
//! use std::io;
//! use std::io::BufReader;
//! use std::str;
//!
//! let mut records = EitherRecords::from(BufReader::new(io::stdin()));
//! while let Some(Ok(record)) = records.next() {
//!     println!("id  : {}", record.id());
//!     println!("desc: {}", record.desc().unwrap_or("none"));
//!     println!("seq : {}", str::from_utf8(record.seq()).unwrap_or("error"));
//!     // Add a default quality in case we have a FASTA record
//!     let default_qual = vec![b'I'; record.seq().len()];
//!     let qual = record.qual().unwrap_or(&default_qual);
//!     println!("qual: {}", str::from_utf8(qual).unwrap_or("error"));
//!     println!("")
//! }
//! ```
//!
//! # Utility Examples
//!
//! ## Type Detection
//!
//! This module provides utilities for detecting if some data is a FASTA/FASTQ.
//!
//! ```
//! use bio::io::fastx::{get_kind, get_kind_file, get_kind_seek, Kind};
//! use std::fs;
//! use std::fs::File;
//! use std::io;
//! use std::io::prelude::*;
//!
//! // From a Read
//!
//! fn from_read() -> io::Result<Kind> {
//!     let reader = io::stdin();
//!     let (mut new_reader, kind) = get_kind(reader)?;
//!     println!("{}", kind);
//!     // Read from start of your old reader
//!     let mut buf = [0u8; 8];
//!     new_reader.read(&mut buf)?;
//!     Ok(kind)
//! }
//!
//! // From a Read + Seek
//!
//! fn from_read_seek() -> io::Result<Kind> {
//!     let mut read_seek = io::Cursor::new(
//!         b">id desc
//! ACTG
//! ",
//!     );
//!
//!     get_kind_seek(&mut read_seek)
//! }
//!
//! // From a file path
//!
//! fn from_file_path() -> io::Result<Kind> {
//!     get_kind_file("foo.fasta")
//! }
//! ```
use anyhow::Context;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::io::BufReader;
use std::io::SeekFrom;
use std::mem;
use std::path::Path;
use thiserror::Error;

use crate::io::{fasta, fastq};
use crate::utils::TextSlice;

macro_rules! passthrough {
    ($name:ident, $t:ty) => {
        fn $name(&self) -> $t {
            self.$name()
        }
    };
}

macro_rules! matchthrough {
    ($name:ident, $t:ty) => {
        fn $name(&self) -> $t {
            match self {
                EitherRecord::FASTA(f) => Record::$name(f),
                EitherRecord::FASTQ(f) => Record::$name(f),
            }
        }
    };
}

pub trait Record {
    fn is_empty(&self) -> bool;
    fn check(&self) -> Result<(), &str>;
    fn id(&self) -> &str;
    fn desc(&self) -> Option<&str>;
    fn seq(&self) -> TextSlice<'_>;
    fn qual(&self) -> Option<&[u8]>;
    fn kind(&self) -> Kind;
}

impl Record for super::fasta::Record {
    passthrough!(is_empty, bool);
    passthrough!(check, Result<(), &str>);
    passthrough!(id, &str);
    passthrough!(desc, Option<&str>);
    passthrough!(seq, TextSlice<'_>);

    fn qual(&self) -> Option<&[u8]> {
        None
    }

    fn kind(&self) -> Kind {
        Kind::FASTA
    }
}

impl Record for super::fastq::Record {
    passthrough!(is_empty, bool);
    passthrough!(check, Result<(), &str>);
    passthrough!(id, &str);
    passthrough!(desc, Option<&str>);
    passthrough!(seq, TextSlice<'_>);

    fn qual(&self) -> Option<&[u8]> {
        Some(self.qual())
    }

    fn kind(&self) -> Kind {
        Kind::FASTQ
    }
}

#[derive(Clone, Display, Debug, Serialize, Deserialize)]
pub enum EitherRecord {
    FASTA(fasta::Record),
    FASTQ(fastq::Record),
}

impl EitherRecord {
    pub fn to_fasta(self) -> fasta::Record {
        return self.into();
    }

    pub fn to_fastq(self, default_qual: u8) -> fastq::Record {
        match self {
            EitherRecord::FASTQ(f) => f,
            EitherRecord::FASTA(f) => {
                let qual = &vec![default_qual; f.seq().len()];
                fastq::Record::with_attrs(f.id(), f.desc(), f.seq(), qual)
            }
        }
    }
}

impl Into<fasta::Record> for EitherRecord {
    fn into(self) -> fasta::Record {
        match self {
            EitherRecord::FASTA(f) => f,
            EitherRecord::FASTQ(f) => fasta::Record::with_attrs(f.id(), f.desc(), f.seq()),
        }
    }
}

impl From<fasta::Record> for EitherRecord {
    fn from(record: fasta::Record) -> Self {
        EitherRecord::FASTA(record)
    }
}

impl From<fastq::Record> for EitherRecord {
    fn from(record: fastq::Record) -> Self {
        EitherRecord::FASTQ(record)
    }
}

impl Record for EitherRecord {
    matchthrough!(is_empty, bool);
    matchthrough!(check, Result<(), &str>);
    matchthrough!(id, &str);
    matchthrough!(desc, Option<&str>);
    matchthrough!(seq, TextSlice<'_>);

    fn qual(&self) -> Option<&[u8]> {
        match &self {
            EitherRecord::FASTA(f) => Record::qual(f),
            EitherRecord::FASTQ(f) => Record::qual(f),
        }
    }

    matchthrough!(kind, Kind);
}

pub trait Records<R: Record, E>: Iterator<Item = Result<R, E>> {}

impl<T: BufRead> Records<fasta::Record, io::Error> for fasta::Records<T> {}

impl<T: BufRead> Records<fastq::Record, fastq::Error> for fastq::Records<T> {}

impl<T: BufRead> Records<EitherRecord, Error> for EitherRecords<T> {}

#[derive(Debug)]
enum EitherRecordsInner<R: BufRead> {
    FASTA(fasta::Records<R>),
    FASTQ(fastq::Records<R>),
}

#[derive(Debug)]
pub struct EitherRecords<R: BufRead> {
    records: Option<EitherRecordsInner<BufReader<io::Chain<io::Cursor<[u8; 1]>, R>>>>,
    reader: Option<R>,
}

impl EitherRecords<BufReader<fs::File>> {
    /// Read from a given file.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(path.as_ref())
            .map(BufReader::new)
            .map(EitherRecords::from)
            .with_context(|| format!("Failed to read fastq from {:#?}", path))
    }
}

impl<R: BufRead> EitherRecords<R> {
    pub fn new(reader: R) -> Self {
        EitherRecords::from(reader)
    }

    pub fn kind(&mut self) -> io::Result<Kind> {
        self.initialize()?;
        match self.records {
            Some(EitherRecordsInner::FASTA(_)) => Ok(Kind::FASTA),
            Some(EitherRecordsInner::FASTQ(_)) => Ok(Kind::FASTQ),
            None => Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                "Data is empty",
            )),
        }
    }

    fn initialize(&mut self) -> io::Result<()> {
        if let Some(reader) = mem::replace(&mut self.reader, None) {
            match get_kind(reader) {
                Err(err) if err.kind() == io::ErrorKind::UnexpectedEof => (),
                Err(err) => return Err(err),
                Ok((reader, Kind::FASTA)) => {
                    self.records = Some(EitherRecordsInner::FASTA(
                        fasta::Reader::new(reader).records(),
                    ))
                }
                Ok((reader, Kind::FASTQ)) => {
                    self.records = Some(EitherRecordsInner::FASTQ(
                        fastq::Reader::new(reader).records(),
                    ))
                }
            }
        }
        Ok(())
    }
}

impl<R: BufRead> Iterator for EitherRecords<R> {
    type Item = Result<EitherRecord>;
    fn next(&mut self) -> Option<Self::Item> {
        if let Err(e) = self.initialize() {
            return Some(Err(Error::IO(e)));
        }
        match &mut self.records {
            Some(EitherRecordsInner::FASTA(r)) => r
                .next()
                .map(|record_res| record_res.map(EitherRecord::FASTA).map_err(Error::IO)),
            Some(EitherRecordsInner::FASTQ(r)) => r
                .next()
                .map(|record_res| record_res.map(EitherRecord::FASTQ).map_err(Error::FASTQ)),
            None => None,
        }
    }
}

impl<R: BufRead> From<R> for EitherRecords<R> {
    fn from(reader: R) -> Self {
        EitherRecords {
            records: None,
            reader: Some(reader),
        }
    }
}

#[derive(Display, Debug, Error)]
pub enum Error {
    IO(io::Error),
    FASTQ(fastq::Error),
}

impl From<io::Error> for Error {
    fn from(err: io::Error) -> Self {
        Error::IO(err)
    }
}

impl From<fastq::Error> for Error {
    fn from(err: fastq::Error) -> Self {
        Error::FASTQ(err)
    }
}

type Result<T, E = Error> = std::result::Result<T, E>;

#[derive(Debug, Eq, PartialEq)]
pub enum Kind {
    FASTQ,
    FASTA,
}

/// Determine whether a [`Read`](Read) is a FASTA or FASTQ.
///
/// This function is a wrapper around [`get_kind_detailed`](get_kind_detailed) for more convenient
/// use if you don't want to:
///
/// - Access your input [`Read`](Read) in the event of an error
/// - Explicitly handle errors resulting from reading from your [`Read`](Read) differently
///   from errors resulting from an invalid FASTA/FASTQ
///
/// You should also only use this function if your [`Read`](Read) does not implement [`Seek`](Seek)
/// othwerwise [`get_kind_seek`](get_kind_seek) is more convenient.
///
/// This method takes ownership of the [`Read`](Read) and returns a new [`Read`](Read)
/// that is equivalent in content to the input [`Read`](Read). This allows you to pass the returned
/// [`Read`]#Read) directly into [`fasta::Reader::new`](fasta::Reader::new) or
/// [`fastq::Reader::new`](fastq::Reader::new).
///
/// # Example
///
/// ```rust
/// use bio::io::fastx::{get_kind, Kind};
/// use bio::io::{fasta, fastq};
/// use std::io;
///
/// fn count_records() -> io::Result<usize> {
///     let (reader, kind) = get_kind(io::stdin())?;
///     match kind {
///         Kind::FASTA => Ok(fasta::Reader::new(reader).records().count()),
///         Kind::FASTQ => Ok(fastq::Reader::new(reader).records().count()),
///     }
/// }
/// ```
pub fn get_kind<R: Read>(reader: R) -> io::Result<(io::Chain<io::Cursor<[u8; 1]>, R>, Kind)> {
    match get_kind_detailed(reader) {
        Ok((reader, Ok(kind))) => Ok((reader, kind)),
        Ok((_, Err(err))) => Err(err),
        Err((_, err)) => Err(err),
    }
}

/// Determine whether a [`Read`](Read) is a FASTA or FASTQ.
///
/// You should only use this function if you would like to:
///
/// - Access your input [`Read`](Read) in the event of an error
/// - Explicitly handle errors resulting from reading from your [`Read`](Read) differently
///   from errors resulting from an invalid FASTA/FASTQ
///
/// Otherwise [`get_kind`](get_kind) will be more convenient.
///
/// You should also only use this function if your [`Read`](Read) does not implement [`Seek`](Seek)
/// othwerwise [`get_kind_seek`](get_kind_seek) is more convenient.
///
/// This method takes ownership of the input [`Read`](Read). It reads from the [`Read`](Read)
/// to determine whether the data is in the FASTA or FASTQ format. If this read fails the
/// function returns an [`Err`](Err) containing a tuple of the input read and the error it
/// encountered (`return Err((read, err))`). Note that in this case there are no guarantee
/// about the position in the returned reader and some data may become unreadable if the position
/// was advanced. Based on the contents of the error the caller can determine how best to recover.
/// If the read is successful this function will always return [`Ok`](Ok) for the outer
/// [`Result`](Result) with a tuple containing a new [`Read`](Read) and a [`Result`](io::Result)
/// containing the [`fastx::Kind`](Kind) if the data was a valid FASTA or FASTQ file or an
/// [`io::Error`](io::Error) if the format was invalid (ex. `Ok((new_reader), Ok(Kind::FASTA))`).
/// The new [`Read`](Read) has the character this function read prepended to it so the
/// resulting [`Read`](Read) should be equivalent in content to the input [`Read`](Read).
///
/// # Example
///
/// ```rust
/// use bio::io::fastx::{get_kind_detailed, Kind};
/// use bio::io::{fasta, fastq};
/// use std::io;
///
/// fn print_type() {
///     let res = get_kind_detailed(io::stdin());
///     match res {
///         Ok((_reader, Ok(Kind::FASTA))) => println!("{}", Kind::FASTA),
///         Ok((_reader, Ok(Kind::FASTQ))) => println!("{}", Kind::FASTQ),
///         Ok((mut reader, Err(e))) => {
///             println!("Error determining FASTA/FASTQ: {}", e);
///             println!("Data:");
///             io::copy(&mut reader, &mut io::stdout());
///         }
///         Err((mut reader, e)) => {
///             println!("Encountered an error while determining FASTA/FASTQ: {}", e);
///             println!("Remaing data:");
///             io::copy(&mut reader, &mut io::stdout());
///         }
///     }
/// }
/// ```
pub fn get_kind_detailed<R: Read>(
    mut reader: R,
) -> std::result::Result<(io::Chain<io::Cursor<[u8; 1]>, R>, io::Result<Kind>), (R, io::Error)> {
    let mut buf = [0];
    if let Err(err) = reader.read_exact(&mut buf) {
        return Err((reader, err));
    }
    let first = char::from(buf[0]);
    let new_reader = io::Cursor::new(buf).chain(reader);

    match first {
        '>' => Ok((new_reader, Ok(Kind::FASTA))),
        '@' => Ok((new_reader, Ok(Kind::FASTQ))),
        _ => Ok((
            new_reader,
            Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "Not a valid FASTA/FASTQ, illegal start character '{}'",
                    first
                ),
            )),
        )),
    }
}

/// Determine whether a [`Read`](Read) + [`Seek`](Seek) is a FASTA or FASTQ.
///
/// The benefit of this this function compared to [`get_kind`](get_kind) is that
/// this function does not take ownership of the [`Read`](Read) so it can
/// be slightly more convenient to use.
pub fn get_kind_seek<R: Read + io::Seek>(reader: &mut R) -> io::Result<Kind> {
    let mut buf = [0];
    reader.read_exact(&mut buf)?;
    reader.seek(SeekFrom::Current(-1))?;
    let first = char::from(buf[0]);

    match first {
        '>' => Ok(Kind::FASTA),
        '@' => Ok(Kind::FASTQ),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "Not a valid FASTA/FASTQ, illegal start character '{}'",
                first
            ),
        )),
    }
}

/// Determine whether a file is a FASTA or FASTQ.
pub fn get_kind_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> io::Result<Kind> {
    fs::File::open(&path).and_then(|mut f| get_kind_seek(&mut f))
}

impl std::fmt::Display for Kind {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Kind::FASTA => write!(f, "{}", "FASTA"),
            Kind::FASTQ => write!(f, "{}", "FASTQ"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use std::io::Cursor;

    const FASTA_FILE: &[u8] = b">id desc
ACCGTAGGCTGA
CCGTAGGCTGAA
CGTAGGCTGAAA
GTAGGCTGAAAA
CCCC
>id2
ATTGTTGTTTTA
ATTGTTGTTTTA
ATTGTTGTTTTA
GGGG
";

    const FASTQ_FILE: &[u8] = b"@id desc
ACCGTAGGCTGA
+
IIIIIIJJJJJJ
";

    const INCOMPLETE_FASTQ_FILE: &[u8] = b"@id desc
ACCGTAGGCTGA
+
";

    struct MockBufReader<R: BufRead> {
        bytes_read: usize,
        error: Option<io::Error>,
        error_at_byte: usize,
        reader: R,
    }

    impl<R: BufRead> MockBufReader<R> {
        fn new(error: io::Error, error_at_byte: usize, reader: R) -> Self {
            MockBufReader {
                bytes_read: 0,
                error: Some(error),
                error_at_byte,
                reader,
            }
        }
    }

    impl<R: BufRead> Read for MockBufReader<R> {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if buf.len() + self.bytes_read >= self.error_at_byte {
                if let Some(err) = mem::replace(&mut self.error, None) {
                    return Err(err);
                }
            }
            let bytes_read = self.reader.read(buf)?;
            self.bytes_read += bytes_read;
            Ok(bytes_read)
        }
    }

    impl<R: BufRead> BufRead for MockBufReader<R> {
        fn fill_buf(&mut self) -> io::Result<&[u8]> {
            Ok(&[])
        }

        fn consume(&mut self, _amt: usize) {}
    }

    #[test]
    fn test_fasta_either_record() {
        let record = EitherRecord::FASTA(fasta::Record::with_attrs("id", Some("desc"), b"ACTG"));
        assert!(matches!(record.kind(), Kind::FASTA));
        assert!(matches!(record.qual(), None));
        let fastq = record.clone().to_fastq(b'I');
        assert_eq!(fastq.id(), "id");
        assert_eq!(fastq.qual(), b"IIII");
        let fasta = record.to_fasta();
        assert_eq!(fasta.id(), "id");
    }

    #[test]
    fn test_fastq_either_record() {
        let record = EitherRecord::FASTQ(fastq::Record::with_attrs(
            "id",
            Some("desc"),
            b"ACTG",
            b"JJJJ",
        ));
        assert!(matches!(record.kind(), Kind::FASTQ));
        assert!(matches!(record.qual(), Some(_)));
        let fastq = record.clone().to_fastq(b'I');
        assert_eq!(fastq.id(), "id");
        assert_eq!(fastq.qual(), b"JJJJ");
        let fasta = record.to_fasta();
        assert_eq!(fasta.id(), "id");
    }

    #[test]
    fn test_records_trait() {
        fn count_records<R: Record, E, I: Records<R, E>>(records: I) -> usize {
            records.count()
        }
        let records = fasta::Reader::new(FASTA_FILE).records();
        let count = count_records(records);
        assert_eq!(count, 2);
    }

    #[test]
    fn test_kind_display() {
        assert_eq!(format!("{}", Kind::FASTA), "FASTA");
        assert_eq!(format!("{}", Kind::FASTQ), "FASTQ");
    }

    #[test]
    fn test_fasta_either_records() {
        let mut records = EitherRecords::new(FASTA_FILE);
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        assert_eq!(records.next().unwrap().unwrap().id(), "id2");
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_fasta_either_records_err() {
        // Error on second read
        // If the error is on the first read it won't test our error conversion as the error
        // will occur during initialization
        let error_at_byte = FASTA_FILE
            .bytes()
            .enumerate()
            .find(move |(i, c)| *i > 0 as usize && c.as_ref().map_or(false, |c| *c == b'>'))
            .unwrap()
            .0
            + 1;
        let reader = MockBufReader::new(
            io::Error::new(io::ErrorKind::Other, "error"),
            error_at_byte,
            FASTA_FILE,
        );
        let mut records = EitherRecords::new(reader);
        assert!(matches!(records.next().unwrap().unwrap_err(), Error::IO(_)));
    }

    #[test]
    fn test_fastq_either_records() {
        let mut records = EitherRecords::new(FASTQ_FILE);
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_fastq_either_records_err() {
        let mut records = EitherRecords::new(INCOMPLETE_FASTQ_FILE);
        assert!(matches!(
            records.next().unwrap().unwrap_err(),
            Error::FASTQ(_)
        ));
    }

    #[test]
    fn test_fasta_either_records_kind() {
        let mut records = EitherRecords::new(FASTA_FILE);
        assert!(matches!(records.kind(), Ok(Kind::FASTA)));
    }

    #[test]
    fn test_fastq_either_records_kind() {
        let mut records = EitherRecords::new(FASTQ_FILE);
        assert!(matches!(records.kind(), Ok(Kind::FASTQ)));
    }

    #[test]
    fn test_empty_either_records_kind() {
        let mut records = EitherRecords::new(b"".as_ref());
        assert!(matches!(records.kind(), Err(_)));
    }

    #[test]
    fn test_empty_either_records() {
        let mut records = EitherRecords::new(b"".as_ref());
        assert!(records.next().is_none());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_invalid_either_records() {
        let mut records = EitherRecords::new(b"(".as_ref());
        assert!(records.next().unwrap().is_err());
        // this second check is intentional
        assert!(records.next().is_none());
    }

    #[test]
    fn test_either_records_from_file() {
        let mut f = fs::File::create("either_records_from_file.fasta").unwrap();
        f.write_all(FASTQ_FILE).unwrap();
        let mut records = EitherRecords::from_file("either_records_from_file.fasta").unwrap();
        assert_eq!(records.next().unwrap().unwrap().id(), "id");
        fs::remove_file("either_records_from_file.fasta").unwrap();
    }

    #[test]
    fn test_get_kind_detailed_read_fasta() {
        let res = get_kind_detailed(FASTA_FILE);
        assert!(matches!(res, Ok((_, Ok(Kind::FASTA)))));
        let mut buf = [0u8; 1];
        res.unwrap().0.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTA_FILE[0]);
    }

    #[test]
    fn test_get_kind_detailed_read_fastq() {
        let (mut reader, kind_res) = get_kind_detailed(FASTQ_FILE).unwrap();
        assert!(matches!(kind_res, Ok(Kind::FASTQ)));
        let mut buf = [0u8; 1];
        reader.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTQ_FILE[0]);
    }

    #[test]
    fn test_get_kind_detailed_read_empty() {
        let (_, err) = get_kind_detailed(Cursor::new(b"")).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_get_kind_detailed_read_invalid() {
        let (mut reader, kind_res) = get_kind_detailed(Cursor::new(b"*")).unwrap();

        assert!(matches!(
            kind_res.unwrap_err().kind(),
            io::ErrorKind::InvalidData
        ));
        let mut buf = [0u8; 1];
        reader.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], b'*');
    }

    #[test]
    fn test_get_kind_detailed_error() {
        let mock_reader =
            MockBufReader::new(io::Error::new(io::ErrorKind::Other, "error"), 0, FASTA_FILE);
        let (mut reader, err) = get_kind_detailed(mock_reader).err().unwrap();
        assert!(matches!(err.kind(), io::ErrorKind::Other));
        let mut buf = [0u8; 1];
        reader.read_exact(&mut buf).unwrap();
        assert_eq!(buf[0], FASTA_FILE[0]);
    }

    #[test]
    fn test_get_kind_seek_fasta() {
        let mut read_seeker = Cursor::new(FASTA_FILE);
        let fastx_kind = get_kind_seek(&mut read_seeker).unwrap();
        assert_eq!(fastx_kind, Kind::FASTA);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_kind_seek_fastq() {
        let mut read_seeker = Cursor::new(FASTQ_FILE);
        let fastq_kind = get_kind_seek(&mut read_seeker).unwrap();
        assert_eq!(fastq_kind, Kind::FASTQ);
        assert_eq!(read_seeker.position(), 0);
    }

    #[test]
    fn test_get_kind_seek_empty() {
        let mut read_seeker = Cursor::new(b"");
        let e = get_kind_seek(&mut read_seeker).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::UnexpectedEof);
    }

    #[test]
    fn test_get_kind_seek_invalid() {
        let mut read_seeker = Cursor::new(b"*");
        let e = get_kind_seek(&mut read_seeker).unwrap_err();
        assert_eq!(e.kind(), io::ErrorKind::InvalidData);
    }

    #[test]
    fn test_get_kind_file() {
        let mut f = fs::File::create("get_kind_file.fasta").unwrap();
        f.write_all(FASTQ_FILE).unwrap();
        let res = get_kind_file("get_kind_file.fasta").unwrap();
        assert_eq!(res, Kind::FASTQ);
        fs::remove_file("get_kind_file.fasta").unwrap();
    }

    #[test]
    fn test_either_record_from_records() {
        let from_fasta = EitherRecord::from(fasta::Record::with_attrs("asd", None, &[]));
        assert_eq!(from_fasta.id(), "asd");
        let from_fastq = EitherRecord::from(fastq::Record::with_attrs("asd", None, &[], &[]));
        assert_eq!(from_fastq.id(), "asd");
    }
}

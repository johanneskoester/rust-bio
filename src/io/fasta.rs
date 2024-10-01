// Copyright 2014-2018 Johannes Köster, Christopher Schröder, Henning Timm.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Structs and trait to read and write files in FASTA format.
//!
//! # Example
//!
//! ## Read
//!
//! In this example, we parse a fasta file from stdin and compute some statistics
//!
//! ```no_run
//! use bio::io::fasta;
//! use std::io;
//!
//! let mut reader = fasta::Reader::new(io::stdin());
//!
//! let mut nb_reads = 0;
//! let mut nb_bases = 0;
//!
//! for result in reader.records() {
//!     let record = result.expect("Error during fasta record parsing");
//!     println!("{}", record.id());
//!
//!     nb_reads += 1;
//!     nb_bases += record.seq().len();
//! }
//!
//! println!("Number of reads: {}", nb_reads);
//! println!("Number of bases: {}", nb_bases);
//! ```
//!
//! We can also use a `while` loop to iterate over records.
//! This is slightly faster than the `for` loop.
//! ```no_run
//! use bio::io::fasta;
//! use std::io;
//! let mut records = fasta::Reader::new(io::stdin()).records();
//!
//! let mut nb_reads = 0;
//! let mut nb_bases = 0;
//!
//! while let Some(Ok(record)) = records.next() {
//!     nb_reads += 1;
//!     nb_bases += record.seq().len();
//! }
//!
//! println!("Number of reads: {}", nb_reads);
//! println!("Number of bases: {}", nb_bases);
//! ```
//!
//! ## Write
//!
//! In this example we generate 10 random sequences with length 100 and write them to stdout.
//!
//! ```
//! use std::io;
//! use bio::io::fasta;
//!
//! let mut seed = 42;
//!
//! let nucleotides = [b'A', b'C', b'G', b'T'];
//!
//! let mut writer = fasta::Writer::new(io::stdout());
//!
//! for _ in 0..10 {
//!     let seq = (0..100).map(|_| {
//!         seed = ((seed ^ seed << 13) ^ seed >> 7) ^ seed << 17; // don't use this random generator
//!         nucleotides[seed % 4]
//!     }).collect::<Vec<u8>>();
//!
//!    writer.write("random", None, seq.as_slice()).expect("Error writing record.");
//! }
//! ```
//!
//! ## Read and Write
//!
//! In this example we filter reads from stdin on sequence length and write them to stdout
//!
//! ```no_run
//! use bio::io::fasta;
//! use bio::io::fasta::FastaRead;
//! use std::io;
//!
//! let mut reader = fasta::Reader::new(io::stdin());
//! let mut writer = fasta::Writer::new(io::stdout());
//! let mut record = fasta::Record::new();
//!
//! while let Ok(()) = reader.read(&mut record) {
//!     if record.is_empty() {
//!         break;
//!     }
//!
//!     if record.seq().len() > 100 {
//!         writer
//!             .write_record(&record)
//!             .ok()
//!             .expect("Error writing record.");
//!     }
//! }
//! ```
//!
//! ## Index
//!
//! Random access to FASTA files is facilitated by [`Index`] and [`IndexedReader`]. The FASTA files
//! must already be indexed with [`samtools faidx`](https://www.htslib.org/doc/faidx.html).
//!
//! In this example, we read in the first 10 bases of the sequence named "chr1".
//!
//! ```rust
//! use bio::io::fasta::IndexedReader;
//! // create dummy files
//! const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
//! const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
//!
//! let seq_name = "chr1";
//! let start: u64 = 0; // start is 0-based, inclusive
//! let stop: u64 = 10; // stop is 0-based, exclusive
//!                     // load the index
//! let mut faidx = IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
//! // move the pointer in the index to the desired sequence and interval
//! faidx
//!     .fetch(seq_name, start, stop)
//!     .expect("Couldn't fetch interval");
//! // read the subsequence defined by the interval into a vector
//! let mut seq = Vec::new();
//! faidx.read(&mut seq).expect("Couldn't read the interval");
//! assert_eq!(seq, b"GTAGGCTGAA");
//! ```

use std::cmp::min;
use std::collections;
use std::convert::AsRef;
use std::fs;
use std::io;
use std::io::prelude::*;
use std::path::Path;

use crate::utils::{Text, TextSlice};
use anyhow::Context;
use std::fmt;

/// Maximum size of temporary buffer used for reading indexed FASTA files.
const MAX_FASTA_BUFFER_SIZE: usize = 512;

/// Trait for FASTA readers.
pub trait FastaRead {
    fn read(&mut self, record: &mut Record) -> io::Result<()>;
}

/// A FASTA reader.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Reader<B> {
    reader: B,
    line: String,
}

impl Reader<io::BufReader<fs::File>> {
    /// Read FASTA from given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: P) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(Reader::new)
            .with_context(|| format!("Failed to read fasta from {:#?}", path))
    }

    /// Read FASTA from give file path and a capacity
    pub fn from_file_with_capacity<P: AsRef<Path> + std::fmt::Debug>(
        capacity: usize,
        path: P,
    ) -> anyhow::Result<Self> {
        fs::File::open(&path)
            .map(|file| Reader::with_capacity(capacity, file))
            .with_context(|| format!("Failed to read fasta from {:#?}", path))
    }
}

impl<R> Reader<io::BufReader<R>>
where
    R: io::Read,
{
    /// Create a new Fasta reader given an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let reader = Reader::new(fasta_file);
    /// # }
    /// ```
    pub fn new(reader: R) -> Self {
        Reader {
            reader: io::BufReader::new(reader),
            line: String::new(),
        }
    }

    /// Create a new Fasta reader given a capacity and an instance of `io::Read`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let reader = Reader::with_capacity(16384, fasta_file);
    /// # }
    /// ```
    pub fn with_capacity(capacity: usize, reader: R) -> Self {
        Reader {
            reader: io::BufReader::with_capacity(capacity, reader),
            line: String::new(),
        }
    }
}

impl<B> Reader<B>
where
    B: io::BufRead,
{
    /// Create a new Fasta reader with an object that implements `io::BufRead`.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// let buffer = io::BufReader::with_capacity(16384, fasta_file);
    /// let reader = Reader::from_bufread(buffer);
    /// # }
    /// ```
    pub fn from_bufread(bufreader: B) -> Self {
        Reader {
            reader: bufreader,
            line: String::new(),
        }
    }

    /// Return an iterator over the records of this Fasta file.
    ///
    /// # Example
    /// ```rust
    /// # use std::io;
    /// # use bio::io::fasta::Reader;
    /// # use bio::io::fasta::Record;
    /// # fn main() {
    /// # const fasta_file: &'static [u8] = b">id desc
    /// # AAAA
    /// # ";
    /// # let reader = Reader::new(fasta_file);
    /// for record in reader.records() {
    ///     let record = record.unwrap();
    ///     assert_eq!(record.id(), "id");
    ///     assert_eq!(record.desc().unwrap(), "desc");
    ///     assert_eq!(record.seq().to_vec(), b"AAAA");
    /// }
    /// # }
    /// ```
    pub fn records(self) -> Records<B> {
        Records {
            reader: self,
            error_has_occured: false,
        }
    }
}

impl<B> FastaRead for Reader<B>
where
    B: io::BufRead,
{
    /// Read the next FASTA record into the given `Record`.
    /// An empty record indicates that no more records can be read.
    ///
    /// Use this method when you want to read records as fast as
    /// possible because it allows the reuse of a `Record` allocation.
    ///
    /// The [records](Reader::records) iterator provides a more ergonomic
    /// approach to accessing FASTA records.
    ///
    /// # Errors
    ///
    /// This function will return an error if the record is incomplete,
    /// syntax is violated or any form of I/O error is encountered.
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fasta::Record;
    /// use bio::io::fasta::{FastaRead, Reader};
    ///
    /// const fasta_file: &'static [u8] = b">id desc
    /// AAAA
    /// ";
    /// let mut reader = Reader::new(fasta_file);
    /// let mut record = Record::new();
    ///
    /// // Check for errors parsing the record
    /// reader
    ///     .read(&mut record)
    ///     .expect("fasta reader: got an io::Error or could not read_line()");
    ///
    /// assert_eq!(record.id(), "id");
    /// assert_eq!(record.desc().unwrap(), "desc");
    /// assert_eq!(record.seq().to_vec(), b"AAAA");
    /// ```
    fn read(&mut self, record: &mut Record) -> io::Result<()> {
        record.clear();
        if self.line.is_empty() {
            self.reader.read_line(&mut self.line)?;
            if self.line.is_empty() {
                return Ok(());
            }
        }

        if !self.line.starts_with('>') {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Expected > at record start.",
            ));
        }
        let mut header_fields = self.line[1..].trim_end().splitn(2, char::is_whitespace);
        record.id = header_fields.next().map(|s| s.to_owned()).unwrap();
        record.desc = header_fields.next().map(|s| s.to_owned());
        loop {
            self.line.clear();
            self.reader.read_line(&mut self.line)?;
            if self.line.is_empty() || self.line.starts_with('>') {
                break;
            }
            record.seq.push_str(self.line.trim_end());
        }

        Ok(())
    }
}

/// A FASTA index as created by SAMtools (.fai).
#[derive(Default, Clone, Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct Index {
    inner: Vec<IndexRecord>,
    name_to_rid: collections::HashMap<String, usize>,
}

impl Index {
    /// Open a FASTA index from a given `io::Read` instance.
    pub fn new<R: io::Read>(fai: R) -> csv::Result<Self> {
        let mut inner = vec![];
        let mut name_to_rid = collections::HashMap::new();

        let mut fai_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
            .from_reader(fai);
        for (rid, row) in fai_reader.deserialize().enumerate() {
            let record: IndexRecord = row?;
            name_to_rid.insert(record.name.clone(), rid);
            inner.push(record);
        }
        Ok(Index { inner, name_to_rid })
    }

    /// Open a FASTA index from a given file path.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: &P) -> anyhow::Result<Self> {
        fs::File::open(path)
            .map_err(csv::Error::from)
            .and_then(Self::new)
            .with_context(|| format!("Failed to read fasta index from {:#?}", path))
    }

    /// Open a FASTA index given the corresponding FASTA file path.
    /// That is, for ref.fasta we expect ref.fasta.fai.
    pub fn with_fasta_file<P: AsRef<Path>>(fasta_path: &P) -> anyhow::Result<Self> {
        let mut fai_path = fasta_path.as_ref().as_os_str().to_owned();
        fai_path.push(".fai");

        Self::from_file(&fai_path)
    }

    /// Return a vector of sequences described in the index.
    pub fn sequences(&self) -> Vec<Sequence> {
        // sort kv pairs by rid to preserve order
        self.inner
            .iter()
            .map(|record| Sequence {
                name: record.name.clone(),
                len: record.len,
            })
            .collect()
    }
}

/// A FASTA reader with an index as created by SAMtools (.fai).
#[derive(Debug)]
pub struct IndexedReader<R: io::Read + io::Seek> {
    reader: io::BufReader<R>,
    pub index: Index,
    fetched_idx: Option<IndexRecord>,
    start: Option<u64>,
    stop: Option<u64>,
}

impl IndexedReader<fs::File> {
    /// Read from a given file path. This assumes the index ref.fasta.fai to be
    /// present for FASTA ref.fasta.
    pub fn from_file<P: AsRef<Path> + std::fmt::Debug>(path: &P) -> anyhow::Result<Self> {
        let index = Index::with_fasta_file(path)?;
        fs::File::open(path)
            .map(|f| Self::with_index(f, index))
            .map_err(csv::Error::from)
            .with_context(|| format!("Failed to read fasta from {:#?}", path))
    }
}

impl<R: io::Read + io::Seek> IndexedReader<R> {
    /// Read from a FASTA and its index, both given as `io::Read`. FASTA has to
    /// be `io::Seek` in addition.
    pub fn new<I: io::Read>(fasta: R, fai: I) -> csv::Result<Self> {
        let index = Index::new(fai)?;
        Ok(IndexedReader {
            reader: io::BufReader::new(fasta),
            index,
            fetched_idx: None,
            start: None,
            stop: None,
        })
    }

    /// Read from a FASTA and its index, the first given as `io::Read`, the
    /// second given as index object.
    pub fn with_index(fasta: R, index: Index) -> Self {
        IndexedReader {
            reader: io::BufReader::new(fasta),
            index,
            fetched_idx: None,
            start: None,
            stop: None,
        }
    }

    /// Fetch an interval from the sequence with the given name for reading.
    ///
    /// `start` and `stop` are 0-based and `stop` is exclusive - i.e. `[start, stop)`
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fasta::IndexedReader;
    /// // create dummy files
    /// const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
    /// const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
    ///
    /// let seq_name = "chr1";
    /// let start: u64 = 0; // start is 0-based, inclusive
    /// let stop: u64 = 10; // stop is 0-based, exclusive
    ///                     // load the index
    /// let mut faidx = IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
    /// // move the pointer in the index to the desired sequence and interval
    /// faidx
    ///     .fetch(seq_name, start, stop)
    ///     .expect("Couldn't fetch interval");
    /// // read the subsequence defined by the interval into a vector
    /// let mut seq = Vec::new();
    /// faidx.read(&mut seq).expect("Couldn't read the interval");
    /// assert_eq!(seq, b"GTAGGCTGAA");
    /// ```
    ///
    /// # Errors
    /// If the `seq_name` does not exist within the index.
    pub fn fetch(&mut self, seq_name: &str, start: u64, stop: u64) -> io::Result<()> {
        let idx = self.idx(seq_name)?;
        self.start = Some(start);
        self.stop = Some(stop);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Fetch an interval from the sequence with the given record index for reading.
    ///
    /// `start` and `stop` are 0-based and `stop` is exclusive - i.e. `[start, stop)`
    ///
    /// # Example
    ///
    /// ```rust
    /// use bio::io::fasta::IndexedReader;
    /// // create dummy files
    /// const FASTA_FILE: &[u8] = b">chr1\nGTAGGCTGAAAA\nCCCC";
    /// const FAI_FILE: &[u8] = b"chr1\t16\t6\t12\t13";
    ///
    /// let rid: usize = 0;
    /// let start: u64 = 0; // start is 0-based, inclusive
    /// let stop: u64 = 10; // stop is 0-based, exclusive
    ///                     // load the index
    /// let mut faidx = IndexedReader::new(std::io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
    /// // move the pointer in the index to the desired sequence and interval
    /// faidx
    ///     .fetch_by_rid(rid, start, stop)
    ///     .expect("Couldn't fetch interval");
    /// // read the subsequence defined by the interval into a vector
    /// let mut seq = Vec::new();
    /// faidx.read(&mut seq).expect("Couldn't read the interval");
    /// assert_eq!(seq, b"GTAGGCTGAA");
    /// ```
    ///
    /// # Errors
    /// If `rid` does not exist within the index.
    pub fn fetch_by_rid(&mut self, rid: usize, start: u64, stop: u64) -> io::Result<()> {
        let idx = self.idx_by_rid(rid)?;
        self.start = Some(start);
        self.stop = Some(stop);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Fetch the whole sequence with the given name for reading.
    pub fn fetch_all(&mut self, seq_name: &str) -> io::Result<()> {
        let idx = self.idx(seq_name)?;
        self.start = Some(0);
        self.stop = Some(idx.len);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Fetch the whole sequence with the given record index for reading.
    pub fn fetch_all_by_rid(&mut self, rid: usize) -> io::Result<()> {
        let idx = self.idx_by_rid(rid)?;
        self.start = Some(0);
        self.stop = Some(idx.len);
        self.fetched_idx = Some(idx);
        Ok(())
    }

    /// Read the fetched sequence into the given vector.
    pub fn read(&mut self, seq: &mut Text) -> io::Result<()> {
        let idx = self.fetched_idx.clone();
        match (idx, self.start, self.stop) {
            (Some(idx), Some(start), Some(stop)) => self.read_into_buffer(idx, start, stop, seq),
            _ => Err(io::Error::new(
                io::ErrorKind::Other,
                "No sequence fetched for reading.",
            )),
        }
    }

    /// Return an iterator yielding the fetched sequence.
    pub fn read_iter(&mut self) -> io::Result<IndexedReaderIterator<'_, R>> {
        let idx = self.fetched_idx.clone();
        match (idx, self.start, self.stop) {
            (Some(idx), Some(start), Some(stop)) => self.read_into_iter(idx, start, stop),
            _ => Err(io::Error::new(
                io::ErrorKind::Other,
                "No sequence fetched for reading.",
            )),
        }
    }

    fn read_into_buffer(
        &mut self,
        idx: IndexRecord,
        start: u64,
        stop: u64,
        seq: &mut Text,
    ) -> io::Result<()> {
        if stop > idx.len {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "FASTA read interval was out of bounds",
            ));
        } else if start > stop {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid query interval",
            ));
        }

        let mut bases_left = stop - start;
        let mut line_offset = self.seek_to(&idx, start)?;

        seq.clear();
        while bases_left > 0 {
            bases_left -= self.read_line(&idx, &mut line_offset, bases_left, seq)?;
        }

        Ok(())
    }

    fn read_into_iter(
        &mut self,
        idx: IndexRecord,
        start: u64,
        stop: u64,
    ) -> io::Result<IndexedReaderIterator<'_, R>> {
        if stop > idx.len {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "FASTA read interval was out of bounds",
            ));
        } else if start > stop {
            return Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid query interval",
            ));
        }

        let bases_left = stop - start;
        let line_offset = self.seek_to(&idx, start)?;
        let capacity = min(
            MAX_FASTA_BUFFER_SIZE,
            min(bases_left, idx.line_bases) as usize,
        );

        Ok(IndexedReaderIterator {
            reader: self,
            record: idx,
            bases_left,
            line_offset,
            buf: Vec::with_capacity(capacity),
            buf_idx: 0,
        })
    }

    /// Return the IndexRecord for the given sequence name or io::Result::Err
    fn idx(&self, seqname: &str) -> io::Result<IndexRecord> {
        match self.index.name_to_rid.get(seqname) {
            Some(rid) => self.idx_by_rid(*rid),
            None => Err(io::Error::new(
                io::ErrorKind::Other,
                format!("Unknown sequence name: {}.", seqname),
            )),
        }
    }

    /// Return the IndexRecord for the given record index or io::Result::Err
    fn idx_by_rid(&self, rid: usize) -> io::Result<IndexRecord> {
        match self.index.inner.get(rid) {
            Some(record) => Ok(record.clone()),
            None => Err(io::Error::new(
                io::ErrorKind::Other,
                "Invalid record index in fasta file.",
            )),
        }
    }

    /// Seek to the given position in the specified FASTA record. The position
    /// of the cursor on the line that the seek ended on is returned.
    fn seek_to(&mut self, idx: &IndexRecord, start: u64) -> io::Result<u64> {
        assert!(start <= idx.len);

        let line_offset = start % idx.line_bases;
        let line_start = start / idx.line_bases * idx.line_bytes;
        let offset = idx.offset + line_start + line_offset;
        self.reader.seek(io::SeekFrom::Start(offset))?;

        Ok(line_offset)
    }

    /// Tries to read up to `bases_left` bases from the current line into `buf`,
    /// returning the actual number of bases read. Depending on the amount of
    /// whitespace per line, the current `line_offset`, and the amount of bytes
    /// returned from `BufReader::fill_buf`, this function may return Ok(0)
    /// multiple times in a row.
    fn read_line(
        &mut self,
        idx: &IndexRecord,
        line_offset: &mut u64,
        bases_left: u64,
        buf: &mut Vec<u8>,
    ) -> io::Result<u64> {
        let (bytes_to_read, bytes_to_keep) = {
            let src = self.reader.fill_buf()?;
            if src.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::UnexpectedEof,
                    "FASTA file is truncated.",
                ));
            }

            let bases_on_line = idx.line_bases - min(idx.line_bases, *line_offset);
            let bases_in_buffer = min(src.len() as u64, bases_on_line);

            let (bytes_to_read, bytes_to_keep) = if bases_in_buffer <= bases_left {
                let bytes_to_read = min(src.len() as u64, idx.line_bytes - *line_offset);

                (bytes_to_read, bases_in_buffer)
            } else {
                (bases_left, bases_left)
            };

            buf.extend_from_slice(&src[..bytes_to_keep as usize]);
            (bytes_to_read, bytes_to_keep)
        };

        self.reader.consume(bytes_to_read as usize);

        assert!(bytes_to_read > 0);
        *line_offset += bytes_to_read;
        if *line_offset >= idx.line_bytes {
            *line_offset = 0;
        }

        Ok(bytes_to_keep)
    }
}

/// Record of a FASTA index.
#[derive(Clone, Eq, PartialEq, Debug, Serialize, Deserialize)]
struct IndexRecord {
    name: String,
    len: u64,
    offset: u64,
    line_bases: u64,
    line_bytes: u64,
}

/// A sequence record returned by the FASTA index.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Sequence {
    pub name: String,
    pub len: u64,
}

#[derive(Debug)]
pub struct IndexedReaderIterator<'a, R: io::Read + io::Seek> {
    reader: &'a mut IndexedReader<R>,
    record: IndexRecord,
    bases_left: u64,
    line_offset: u64,
    buf: Vec<u8>,
    buf_idx: usize,
}

impl<'a, R: io::Read + io::Seek + 'a> IndexedReaderIterator<'a, R> {
    fn fill_buffer(&mut self) -> io::Result<()> {
        assert!(self.bases_left > 0);

        self.buf.clear();
        let bases_to_read = min(self.buf.capacity() as u64, self.bases_left);

        // May loop one or more times; see IndexedReader::read_line.
        while self.buf.is_empty() {
            self.bases_left -= self.reader.read_line(
                &self.record,
                &mut self.line_offset,
                bases_to_read,
                &mut self.buf,
            )?;
        }

        self.buf_idx = 0;
        Ok(())
    }
}

impl<'a, R: io::Read + io::Seek + 'a> Iterator for IndexedReaderIterator<'a, R> {
    type Item = io::Result<u8>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buf_idx < self.buf.len() {
            let item = Some(Ok(self.buf[self.buf_idx]));
            self.buf_idx += 1;
            item
        } else if self.bases_left > 0 {
            if let Err(e) = self.fill_buffer() {
                self.bases_left = 0;
                self.buf_idx = self.buf.len();

                return Some(Err(e));
            }

            self.buf_idx = 1;
            Some(Ok(self.buf[0]))
        } else {
            None
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let hint = self.bases_left as usize + (self.buf.len() - self.buf_idx);

        (hint, Some(hint))
    }
}

/// A Fasta writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    writer: io::BufWriter<W>,
    linewrap: Option<usize>,
}

impl Writer<fs::File> {
    /// Write to the given file path.
    #[allow(clippy::wrong_self_convention)]
    pub fn to_file<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        fs::File::create(path).map(Writer::new)
    }

    /// Write to the given file path and a buffer capacity
    pub fn to_file_with_capacity<P: AsRef<Path>>(capacity: usize, path: P) -> io::Result<Self> {
        fs::File::create(path).map(|file| Writer::with_capacity(capacity, file))
    }
}

impl<W: io::Write> Writer<W> {
    /// Create a new Fasta writer.
    pub fn new(writer: W) -> Self {
        Writer {
            writer: io::BufWriter::new(writer),
            linewrap: None,
        }
    }

    /// Create a new Fasta writer with a capacity of write buffer
    pub fn with_capacity(capacity: usize, writer: W) -> Self {
        Writer {
            writer: io::BufWriter::with_capacity(capacity, writer),
            linewrap: None,
        }
    }

    /// Create a new Fasta writer with a given BufWriter
    pub fn from_bufwriter(bufwriter: io::BufWriter<W>) -> Self {
        Writer {
            writer: bufwriter,
            linewrap: None,
        }
    }

    /// Directly write a [`fasta::Record`](struct.Record.html).
    ///
    /// # Errors
    /// If there is an issue writing to the `Writer`.
    ///
    /// # Examples
    /// ```rust
    /// use bio::io::fasta::{Record, Writer};
    /// use std::fs;
    /// use std::io;
    /// use std::path::Path;
    ///
    /// let path = Path::new("test.fa");
    /// let file = fs::File::create(path).unwrap();
    /// {
    ///     let handle = io::BufWriter::new(file);
    ///     let mut writer = Writer::new(handle);
    ///     let record = Record::with_attrs("id", Some("desc"), b"ACGT");
    ///
    ///     let write_result = writer.write_record(&record);
    ///     assert!(write_result.is_ok());
    /// }
    ///
    /// let actual = fs::read_to_string(path).unwrap();
    /// let expected = ">id desc\nACGT\n";
    ///
    /// assert!(fs::remove_file(path).is_ok());
    /// assert_eq!(actual, expected)
    /// ```

    pub fn write_record(&mut self, record: &Record) -> io::Result<()> {
        self.write(record.id(), record.desc(), record.seq())
    }

    /// Directly write a [`fasta::Record`](struct.Record.html) width a fixed line length.
    ///
    /// # Errors
    /// If there is an issue writing to the `Writer`.
    ///
    /// # Examples
    /// ```rust
    /// use bio::io::fasta::{Record, Writer};
    /// use std::fs;
    /// use std::io;
    /// use std::path::Path;
    ///
    /// let path = Path::new("test.fa");
    /// let file = fs::File::create(path).unwrap();
    /// {
    ///     let handle = io::BufWriter::new(file);
    ///     let mut writer = Writer::new(handle);
    ///     
    ///     // For demonstration width is 4 chars, use 50, 60 or 70 instead for production
    ///     writer.set_linewrap(Some(4));
    ///
    ///     let record = Record::with_attrs("id", Some("desc"), b"ACGTACGT");
    ///     let write_result = writer.write_record(&record);
    ///     assert!(write_result.is_ok());
    /// }
    ///
    /// let actual = fs::read_to_string(path).unwrap();
    /// let expected = ">id desc\nACGT\nACGT\n";
    ///
    /// assert!(fs::remove_file(path).is_ok());
    /// assert_eq!(actual, expected)
    /// ```

    pub fn set_linewrap(&mut self, linewrap: Option<usize>) {
        self.linewrap = linewrap
    }

    pub fn write_record_header(&mut self, id: &str, desc: Option<&str>) -> io::Result<()> {
        self.writer.write_all(b">")?;
        self.writer.write_all(id.as_bytes())?;
        if let Some(desc) = desc {
            self.writer.write_all(b" ")?;
            self.writer.write_all(desc.as_bytes())?;
        }
        self.writer.write_all(b"\n")?;

        Ok(())
    }

    /// Write a Fasta record with given id, optional description and sequence.
    pub fn write(&mut self, id: &str, desc: Option<&str>, seq: TextSlice<'_>) -> io::Result<()> {
        self.write_record_header(id, desc)?;
        if self.linewrap == None {
            self.writer.write_all(seq)?;
            self.writer.write_all(b"\n")?;
            Ok(())
        } else {
            // Write Fasta lines with a given linewrap instead of in a single line
            seq.chunks(self.linewrap.unwrap())
                .try_for_each(|chunk| -> io::Result<()> {
                    self.writer.write_all(chunk)?;
                    self.writer.write_all(b"\n")?;

                    Ok(())
                })
        }
    }

    /// Flush the writer, ensuring that everything is written.
    pub fn flush(&mut self) -> io::Result<()> {
        self.writer.flush()
    }
}

/// A FASTA record.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Record {
    id: String,
    desc: Option<String>,
    seq: String,
}

impl Record {
    /// Create a new instance.
    pub fn new() -> Self {
        Record {
            id: String::new(),
            desc: None,
            seq: String::new(),
        }
    }

    /// Create a new `Record` from given attributes.
    ///
    /// # Examples
    /// ```rust
    /// use bio::io::fasta::Record;
    ///
    /// let read_id = "read1";
    /// let description = Some("sampleid=foobar");
    /// let sequence = b"ACGT";
    /// let record = Record::with_attrs(read_id, description, sequence);
    ///
    /// assert_eq!(">read1 sampleid=foobar\nACGT\n", record.to_string())
    /// ```
    pub fn with_attrs(id: &str, desc: Option<&str>, seq: TextSlice<'_>) -> Self {
        let desc = desc.map(|desc| desc.to_owned());
        Record {
            id: id.to_owned(),
            desc,
            seq: String::from_utf8(seq.to_vec()).unwrap(),
        }
    }

    /// Check if record is empty.
    pub fn is_empty(&self) -> bool {
        self.id.is_empty() && self.desc.is_none() && self.seq.is_empty()
    }

    /// Check validity of Fasta record.
    pub fn check(&self) -> Result<(), &str> {
        if self.id().is_empty() {
            return Err("Expecting id for Fasta record.");
        }
        if !self.seq.is_ascii() {
            return Err("Non-ascii character found in sequence.");
        }

        Ok(())
    }

    /// Return the id of the record.
    pub fn id(&self) -> &str {
        self.id.as_ref()
    }

    /// Return descriptions if present.
    pub fn desc(&self) -> Option<&str> {
        match self.desc.as_ref() {
            Some(desc) => Some(desc),
            None => None,
        }
    }

    /// Return the sequence of the record.
    pub fn seq(&self) -> TextSlice<'_> {
        self.seq.as_bytes()
    }

    /// Clear the record.
    fn clear(&mut self) {
        self.id.clear();
        self.desc = None;
        self.seq.clear();
    }
}

impl fmt::Display for Record {
    /// Allows for using `Record` in a given formatter `f`. In general this is for
    /// creating a `String` representation of a `Record` and, optionally, writing it to
    /// a file.
    ///
    /// # Errors
    /// Returns [`std::fmt::Error`](https://doc.rust-lang.org/std/fmt/struct.Error.html)
    /// if there is an issue formatting to the stream.
    ///
    /// # Examples
    ///
    /// Read in a Fasta `Record` and create a `String` representation of it.
    ///
    /// ```rust
    /// use bio::io::fasta::Reader;
    /// use std::fmt::Write;
    /// // create a "fake" fasta file
    /// let fasta: &'static [u8] = b">id comment1 comment2\nACGT\n";
    /// let mut records = Reader::new(fasta).records().map(|r| r.unwrap());
    /// let record = records.next().unwrap();
    ///
    /// let mut actual = String::new();
    /// // populate `actual` with a string representation of our record
    /// write!(actual, "{}", record).unwrap();
    ///
    /// let expected = std::str::from_utf8(fasta).unwrap();
    ///
    /// assert_eq!(actual, expected)
    /// ```
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let header = match self.desc() {
            Some(d) => format!("{} {}", self.id().to_owned(), d),
            None => self.id().to_owned(),
        };
        write!(
            f,
            ">{}\n{}\n",
            header,
            std::str::from_utf8(self.seq()).unwrap(),
        )
    }
}

/// An iterator over the records of a Fasta file.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct Records<B>
where
    B: io::BufRead,
{
    reader: Reader<B>,
    error_has_occured: bool,
}

impl<B> Iterator for Records<B>
where
    B: io::BufRead,
{
    type Item = io::Result<Record>;

    fn next(&mut self) -> Option<io::Result<Record>> {
        if self.error_has_occured {
            None
        } else {
            let mut record = Record::new();
            match self.reader.read(&mut record) {
                Ok(()) if record.is_empty() => None,
                Ok(()) => Some(Ok(record)),
                Err(err) => {
                    self.error_has_occured = true;
                    Some(Err(err))
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fmt::Write as FmtWrite;
    use std::io;

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
    const FAI_FILE: &[u8] = b"id\t52\t9\t12\t13
id2\t40\t71\t12\t13
";

    const TRUNCATED_FASTA: &[u8] = b">id desc\nACCGTAGGCTGA";

    const FASTA_FILE_CRLF: &[u8] = b">id desc\r
ACCGTAGGCTGA\r
CCGTAGGCTGAA\r
CGTAGGCTGAAA\r
GTAGGCTGAAAA\r
CCCC\r
>id2\r
ATTGTTGTTTTA\r
ATTGTTGTTTTA\r
ATTGTTGTTTTA\r
GGGG\r
";
    const FAI_FILE_CRLF: &[u8] = b"id\t52\t10\t12\t14\r
id2\t40\t78\t12\t14\r
";

    const FASTA_FILE_NO_TRAILING_LF: &[u8] = b">id desc
GTAGGCTGAAAA
CCCC";
    const FAI_FILE_NO_TRAILING_LF: &[u8] = b"id\t16\t9\t12\t13";

    const WRITE_FASTA_FILE: &[u8] = b">id desc
ACCGTAGGCTGA
>id2
ATTGTTGTTTTA
";
    const WRITE_FASTA_FILE_WIDTH: &[u8] = b">id desc
ACCG
TAGG
CTGA
>id2
ATTG
TTGT
TTTA
";

    struct ReaderMock {
        seek_fails: bool,
        read_fails: bool,
    }

    impl Read for ReaderMock {
        fn read(&mut self, buf: &mut [u8]) -> io::Result<usize> {
            if self.read_fails {
                Err(io::Error::new(io::ErrorKind::Other, "Read set to fail"))
            } else {
                Ok(buf.len())
            }
        }
    }

    impl Seek for ReaderMock {
        fn seek(&mut self, pos: io::SeekFrom) -> io::Result<u64> {
            if let io::SeekFrom::Start(pos) = pos {
                if self.seek_fails {
                    Err(io::Error::new(io::ErrorKind::Other, "Seek set to fail"))
                } else {
                    Ok(pos)
                }
            } else {
                unimplemented!();
            }
        }
    }

    #[test]
    fn test_reader() {
        let reader = Reader::new(FASTA_FILE);
        let ids = ["id", "id2"];
        let descs = [Some("desc"), None];
        let seqs: [&[u8]; 2] = [
            b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC",
            b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG",
        ];

        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
        }

        let reader = Reader::with_capacity(100, FASTA_FILE);

        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
        }

        let reader = Reader::from_bufread(io::BufReader::new(FASTA_FILE));

        for (i, r) in reader.records().enumerate() {
            let record = r.expect("Error reading record");
            assert_eq!(record.check(), Ok(()));
            assert_eq!(record.id(), ids[i]);
            assert_eq!(record.desc(), descs[i]);
            assert_eq!(record.seq(), seqs[i]);
        }
    }

    #[test]
    fn test_faread_trait() {
        let path = "genome.fa.gz";
        let mut fa_reader: Box<dyn FastaRead> = match path.ends_with(".gz") {
            true => Box::new(Reader::new(io::BufReader::new(FASTA_FILE))),
            false => Box::new(Reader::new(FASTA_FILE)),
        };
        // The read method can be called, since it is implemented by
        // FQRead. Right now, the records method would not work.
        let mut record = Record::new();
        fa_reader.read(&mut record).unwrap();
        // Check if the returned result is correct.
        assert_eq!(record.check(), Ok(()));
        assert_eq!(record.id(), "id");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(
            record.seq().to_vec(),
            b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC".to_vec()
        );
    }

    #[test]
    fn test_reader_wrong_header() {
        let mut reader = Reader::new(&b"!test\nACGTA\n"[..]);
        let mut record = Record::new();
        assert!(
            reader.read(&mut record).is_err(),
            "read() should return Err if FASTA header is malformed"
        );
    }

    #[test]
    fn test_reader_no_id() {
        let mut reader = Reader::new(&b">\nACGTA\n"[..]);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        assert!(
            record.check().is_err(),
            "check() should return Err if FASTA header is empty"
        );
    }

    #[test]
    fn test_reader_non_ascii_sequence() {
        let mut reader = Reader::new(&b">id\nACGTA\xE2\x98\xB9AT\n"[..]);
        let mut record = Record::new();
        reader.read(&mut record).unwrap();
        assert!(
            record.check().is_err(),
            "check() should return Err if FASTA sequence is not ASCII"
        );
    }

    #[test]
    fn test_reader_read_fails() {
        let mut reader = Reader::new(ReaderMock {
            seek_fails: false,
            read_fails: true,
        });
        let mut record = Record::new();
        assert!(
            reader.read(&mut record).is_err(),
            "read() should return Err if Read::read fails"
        );
    }

    #[test]
    fn test_reader_read_fails_iter() {
        let reader = Reader::new(ReaderMock {
            seek_fails: false,
            read_fails: true,
        });
        let mut records = reader.records();

        assert!(
            records.next().unwrap().is_err(),
            "next() should return Err if Read::read fails"
        );
        assert!(
            records.next().is_none(),
            "next() should return None after error has occurred"
        );
    }

    #[test]
    fn test_reader_from_file_path_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.fasta");
        let error = Reader::from_file(path)
            .unwrap_err()
            .downcast::<String>()
            .unwrap();

        assert_eq!(&error, "Failed to read fasta from \"/I/dont/exist.fasta\"")
    }

    #[test]
    fn test_record_with_attrs_without_description() {
        let record = Record::with_attrs("id_str", None, b"ATGCGGG");
        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), None);
        assert_eq!(record.seq(), b"ATGCGGG");
    }

    #[test]
    fn test_record_with_attrs_with_description() {
        let record = Record::with_attrs("id_str", Some("desc"), b"ATGCGGG");
        assert_eq!(record.id(), "id_str");
        assert_eq!(record.desc(), Some("desc"));
        assert_eq!(record.seq(), b"ATGCGGG");
    }

    #[test]
    fn test_index_sequences() {
        let reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();

        let sequences = reader.index.sequences();
        assert_eq!(sequences.len(), 2);
        assert_eq!(
            sequences[0],
            Sequence {
                name: "id".into(),
                len: 52,
            }
        );
        assert_eq!(
            sequences[1],
            Sequence {
                name: "id2".into(),
                len: 40,
            }
        );
    }

    #[test]
    fn test_indexed_reader() {
        _test_indexed_reader(FASTA_FILE, FAI_FILE, _read_buffer);
        _test_indexed_reader_truncated(_read_buffer);
        _test_indexed_reader_extreme_whitespace(_read_buffer);
    }

    #[test]
    fn test_indexed_reader_crlf() {
        _test_indexed_reader(FASTA_FILE_CRLF, FAI_FILE_CRLF, _read_buffer);
    }

    #[test]
    fn test_indexed_reader_iter() {
        _test_indexed_reader(FASTA_FILE, FAI_FILE, _read_iter);
        _test_indexed_reader_truncated(_read_iter);
        _test_indexed_reader_extreme_whitespace(_read_iter);
    }

    #[test]
    fn test_indexed_reader_iter_crlf() {
        _test_indexed_reader(FASTA_FILE_CRLF, FAI_FILE_CRLF, _read_iter);
    }

    fn _test_indexed_reader<'a, F>(fasta: &'a [u8], fai: &'a [u8], read: F)
    where
        F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, &str, u64, u64) -> io::Result<Vec<u8>>,
    {
        let mut reader = IndexedReader::new(io::Cursor::new(fasta), fai).unwrap();

        // Test reading various substrings of the sequence
        assert_eq!(read(&mut reader, "id", 1, 5).unwrap(), b"CCGT");
        assert_eq!(
            read(&mut reader, "id", 1, 31).unwrap(),
            b"CCGTAGGCTGACCGTAGGCTGAACGTAGGC"
        );
        assert_eq!(read(&mut reader, "id", 13, 23).unwrap(), b"CGTAGGCTGA");
        assert_eq!(
            read(&mut reader, "id", 36, 52).unwrap(),
            b"GTAGGCTGAAAACCCC"
        );
        assert_eq!(
            read(&mut reader, "id2", 12, 40).unwrap(),
            b"ATTGTTGTTTTAATTGTTGTTTTAGGGG"
        );
        assert_eq!(read(&mut reader, "id2", 12, 12).unwrap(), b"");
        assert_eq!(read(&mut reader, "id2", 12, 13).unwrap(), b"A");
        // Minimal sequence spanning new-line
        assert_eq!(read(&mut reader, "id", 11, 13).unwrap(), b"AC");

        assert!(read(&mut reader, "id2", 12, 11).is_err());
        assert!(read(&mut reader, "id2", 12, 1000).is_err());
        assert!(read(&mut reader, "id3", 0, 1).is_err());
    }

    fn _test_indexed_reader_truncated<'a, F>(read: F)
    where
        F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, &str, u64, u64) -> io::Result<Vec<u8>>,
    {
        let mut reader = IndexedReader::new(io::Cursor::new(TRUNCATED_FASTA), FAI_FILE).unwrap();

        assert_eq!(read(&mut reader, "id", 0, 12).unwrap(), b"ACCGTAGGCTGA");
        assert!(read(&mut reader, "id", 0, 13).is_err()); // read past EOF
        assert!(read(&mut reader, "id", 36, 52).is_err()); // seek and read past EOF
        assert!(read(&mut reader, "id2", 12, 40).is_err()); // seek and read past EOF
    }

    fn _test_indexed_reader_extreme_whitespace<F>(read: F)
    where
        F: Fn(&mut IndexedReader<io::Cursor<Vec<u8>>>, &str, u64, u64) -> io::Result<Vec<u8>>,
    {
        // Test to exercise the case where we cannot consume all whitespace at once. More than
        // DEFAULT_BUF_SIZE (a non-public constant set to 8 * 1024) whitespace is used to ensure
        // that it can't all fit in the BufReader at once.
        let mut seq = Vec::new();
        seq.push(b'A');
        seq.resize(10000, b' ');
        seq.push(b'B');

        let fasta = io::Cursor::new(seq);
        let fai = io::Cursor::new(Vec::from(&b"id\t2\t0\t1\t10000"[..]));
        let mut reader = IndexedReader::new(fasta, fai).unwrap();

        assert_eq!(read(&mut reader, "id", 0, 2).unwrap(), b"AB");
    }

    fn _read_buffer<T>(
        reader: &mut IndexedReader<T>,
        seqname: &str,
        start: u64,
        stop: u64,
    ) -> io::Result<Vec<u8>>
    where
        T: Seek + Read,
    {
        let mut seq = vec![];
        reader.fetch(seqname, start, stop)?;
        reader.read(&mut seq)?;

        Ok(seq)
    }

    fn _read_iter<T>(
        reader: &mut IndexedReader<T>,
        seqname: &str,
        start: u64,
        stop: u64,
    ) -> io::Result<Vec<u8>>
    where
        T: Seek + Read,
    {
        let mut seq = vec![];
        reader.fetch(seqname, start, stop)?;
        for nuc in reader.read_iter()? {
            seq.push(nuc?);
        }

        Ok(seq)
    }

    #[test]
    fn test_indexed_reader_all() {
        _test_indexed_reader_all(FASTA_FILE, FAI_FILE, _read_buffer_all);
    }

    #[test]
    fn test_indexed_reader_crlf_all() {
        _test_indexed_reader_all(FASTA_FILE_CRLF, FAI_FILE_CRLF, _read_buffer_all);
    }

    #[test]
    fn test_indexed_reader_iter_all() {
        _test_indexed_reader_all(FASTA_FILE, FAI_FILE, _read_iter_all);
    }

    #[test]
    fn test_indexed_reader_iter_crlf_all() {
        _test_indexed_reader_all(FASTA_FILE_CRLF, FAI_FILE_CRLF, _read_iter_all);
    }

    fn _test_indexed_reader_all<'a, F>(fasta: &'a [u8], fai: &'a [u8], read: F)
    where
        F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, &str) -> io::Result<Vec<u8>>,
    {
        let mut reader = IndexedReader::new(io::Cursor::new(fasta), fai).unwrap();

        assert_eq!(
            read(&mut reader, "id").unwrap(),
            &b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC"[..]
        );
        assert_eq!(
            read(&mut reader, "id2").unwrap(),
            &b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG"[..]
        );
    }

    fn _read_buffer_all<T>(reader: &mut IndexedReader<T>, seqname: &str) -> io::Result<Vec<u8>>
    where
        T: Seek + Read,
    {
        let mut seq = vec![];
        reader.fetch_all(seqname)?;
        reader.read(&mut seq)?;

        Ok(seq)
    }

    fn _read_iter_all<T>(reader: &mut IndexedReader<T>, seqname: &str) -> io::Result<Vec<u8>>
    where
        T: Seek + Read,
    {
        let mut seq = vec![];
        reader.fetch_all(seqname)?;
        for nuc in reader.read_iter()? {
            seq.push(nuc?);
        }

        Ok(seq)
    }

    #[test]
    fn test_indexed_reader_by_rid_all() {
        _test_indexed_reader_by_rid_all(FASTA_FILE, FAI_FILE, _read_buffer_by_rid_all);
    }

    #[test]
    fn test_indexed_reader_crlf_by_rid_all() {
        _test_indexed_reader_by_rid_all(FASTA_FILE_CRLF, FAI_FILE_CRLF, _read_buffer_by_rid_all);
    }

    #[test]
    fn test_indexed_reader_iter_by_rid_all() {
        _test_indexed_reader_by_rid_all(FASTA_FILE, FAI_FILE, _read_iter_by_rid_all);
    }

    #[test]
    fn test_indexed_reader_iter_crlf_by_rid_all() {
        _test_indexed_reader_by_rid_all(FASTA_FILE_CRLF, FAI_FILE_CRLF, _read_iter_by_rid_all);
    }

    fn _test_indexed_reader_by_rid_all<'a, F>(fasta: &'a [u8], fai: &'a [u8], read: F)
    where
        F: Fn(&mut IndexedReader<io::Cursor<&'a [u8]>>, usize) -> io::Result<Vec<u8>>,
    {
        let mut reader = IndexedReader::new(io::Cursor::new(fasta), fai).unwrap();

        assert_eq!(
            read(&mut reader, 0).unwrap(),
            &b"ACCGTAGGCTGACCGTAGGCTGAACGTAGGCTGAAAGTAGGCTGAAAACCCC"[..]
        );
        assert_eq!(
            read(&mut reader, 1).unwrap(),
            &b"ATTGTTGTTTTAATTGTTGTTTTAATTGTTGTTTTAGGGG"[..]
        );
    }

    fn _read_buffer_by_rid_all<T>(
        reader: &mut IndexedReader<T>,
        seq_index: usize,
    ) -> io::Result<Vec<u8>>
    where
        T: Seek + Read,
    {
        let mut seq = vec![];
        reader.fetch_all_by_rid(seq_index)?;
        reader.read(&mut seq)?;

        Ok(seq)
    }

    fn _read_iter_by_rid_all<T>(
        reader: &mut IndexedReader<T>,
        seq_index: usize,
    ) -> io::Result<Vec<u8>>
    where
        T: Seek + Read,
    {
        let mut seq = vec![];
        reader.fetch_all_by_rid(seq_index)?;
        for nuc in reader.read_iter()? {
            seq.push(nuc?);
        }

        Ok(seq)
    }

    #[test]
    fn test_indexed_reader_iter_size_hint() {
        let mut reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
        reader.fetch("id", 2, 4).unwrap();
        let mut iterator = reader.read_iter().unwrap();

        assert_eq!(iterator.size_hint(), (2, Some(2)));
        assert_eq!(iterator.next().unwrap().unwrap(), b'C');
        assert_eq!(iterator.size_hint(), (1, Some(1)));
        assert_eq!(iterator.next().unwrap().unwrap(), b'G');
        assert_eq!(iterator.size_hint(), (0, Some(0)));
        assert!(iterator.next().is_none());
        assert_eq!(iterator.size_hint(), (0, Some(0)));
    }

    #[test]
    fn test_indexed_reader_reused_buffer() {
        let mut reader = IndexedReader::new(io::Cursor::new(FASTA_FILE), FAI_FILE).unwrap();
        let mut seq = Vec::new();

        reader.fetch("id", 1, 5).unwrap();
        reader.read(&mut seq).unwrap();
        assert_eq!(seq, b"CCGT");

        reader.fetch("id", 13, 23).unwrap();
        reader.read(&mut seq).unwrap();
        assert_eq!(seq, b"CGTAGGCTGA");
    }

    #[test]
    fn test_indexed_reader_no_trailing_lf() {
        let mut reader = IndexedReader::new(
            io::Cursor::new(FASTA_FILE_NO_TRAILING_LF),
            FAI_FILE_NO_TRAILING_LF,
        )
        .unwrap();
        let mut seq = Vec::new();

        reader.fetch("id", 0, 16).unwrap();
        reader.read(&mut seq).unwrap();
        assert_eq!(seq, b"GTAGGCTGAAAACCCC");
    }

    #[test]
    fn test_indexed_reader_bad_reader() {
        let bad_reader = ReaderMock {
            seek_fails: false,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut seq = Vec::new();
        reader.fetch("id", 0, 10).unwrap();
        assert!(reader.read(&mut seq).is_ok())
    }

    #[test]
    fn test_indexed_reader_read_seek_fails() {
        let bad_reader = ReaderMock {
            seek_fails: true,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut seq = Vec::new();
        reader.fetch("id", 0, 10).unwrap();
        assert!(reader.read(&mut seq).is_err());
    }

    #[test]
    fn test_indexed_reader_read_read_fails() {
        let bad_reader = ReaderMock {
            seek_fails: false,
            read_fails: true,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        let mut seq = Vec::new();
        reader.fetch("id", 0, 10).unwrap();
        assert!(reader.read(&mut seq).is_err());
    }

    #[test]
    fn test_indexed_reader_iter_seek_fails() {
        let bad_reader = ReaderMock {
            seek_fails: true,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        reader.fetch("id", 0, 10).unwrap();
        assert!(reader.read_iter().is_err());
    }

    #[test]
    fn test_indexed_reader_iter_read_fails() {
        let bad_reader = ReaderMock {
            seek_fails: false,
            read_fails: true,
        };
        let mut reader = IndexedReader::new(bad_reader, FAI_FILE).unwrap();
        reader.fetch("id", 0, 10).unwrap();
        let mut iterator = reader.read_iter().unwrap();
        assert!(iterator.next().unwrap().is_err());
        assert!(
            iterator.next().is_none(),
            "next() should return none after error has occurred"
        );
    }

    #[test]
    fn test_indexed_reader_no_fetch_read_fails() {
        let reader = ReaderMock {
            seek_fails: false,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(reader, FAI_FILE).unwrap();
        let mut seq = vec![];
        assert!(reader.read(&mut seq).is_err());
    }

    #[test]
    fn test_indexed_reader_no_fetch_read_iter_fails() {
        let reader = ReaderMock {
            seek_fails: false,
            read_fails: false,
        };
        let mut reader = IndexedReader::new(reader, FAI_FILE).unwrap();
        assert!(reader.read_iter().is_err());
    }

    #[test]
    fn test_writer() {
        let mut writer = Writer::new(Vec::new());
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE);

        let mut writer = Writer::with_capacity(100, Vec::new());
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE);

        let mut writer = Writer::from_bufwriter(std::io::BufWriter::with_capacity(100, Vec::new()));
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE);
    }

    #[test]
    fn test_display_record_no_desc_id_without_space_after() {
        let fasta: &'static [u8] = b">id\nACGT\n";
        let mut records = Reader::new(fasta).records().map(|r| r.unwrap());
        let record = records.next().unwrap();
        let mut actual = String::new();
        write!(actual, "{}", record).unwrap();

        let expected = std::str::from_utf8(fasta).unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_display_record_with_desc_id_has_space_between_id_and_desc() {
        let fasta: &'static [u8] = b">id comment1 comment2\nACGT\n";
        let mut records = Reader::new(fasta).records().map(|r| r.unwrap());
        let record = records.next().unwrap();
        let mut actual = String::new();
        write!(actual, "{}", record).unwrap();

        let expected = std::str::from_utf8(fasta).unwrap();

        assert_eq!(actual, expected)
    }

    #[test]
    fn test_index_record_idx_by_rid_invalid_index_returns_error() {
        let reader = ReaderMock {
            seek_fails: false,
            read_fails: false,
        };
        let index_reader = IndexedReader::new(reader, FAI_FILE).unwrap();

        let actual = index_reader.idx_by_rid(99999).unwrap_err();
        let expected = io::Error::new(io::ErrorKind::Other, "Invalid record index in fasta file.");

        assert_eq!(actual.kind(), expected.kind());
        assert_eq!(actual.to_string(), expected.to_string())
    }

    #[test]
    fn test_index_record_fetch_by_rid_second_index_returns_second_record() {
        let reader = ReaderMock {
            seek_fails: false,
            read_fails: false,
        };
        let mut index_reader = IndexedReader::new(reader, FAI_FILE).unwrap();

        let actual = index_reader.fetch_by_rid(1, 1, 3);

        assert!(actual.is_ok());
        assert_eq!(
            index_reader.fetched_idx,
            Some(IndexRecord {
                name: String::from("id2"),
                len: 40,
                offset: 71,
                line_bases: 12,
                line_bytes: 13
            })
        )
    }

    #[test]
    fn test_writer_to_file_dir_doesnt_exist_returns_err() {
        let path = Path::new("/I/dont/exist.fa");

        let actual = Writer::to_file(path).unwrap_err();
        let expected = io::Error::new(io::ErrorKind::NotFound, "foo");

        assert_eq!(actual.kind(), expected.kind());
    }

    #[test]
    fn test_writer_to_file_dir_exists_returns_ok() {
        let file = tempfile::NamedTempFile::new().expect("Could not create temp file");
        let path = file.path();

        assert!(Writer::to_file(path).is_ok());
        assert!(Writer::to_file_with_capacity(100, path).is_ok());
    }

    #[test]
    fn test_write_record() {
        let path = Path::new("test.fa");
        let file = fs::File::create(path).unwrap();
        {
            let handle = io::BufWriter::new(file);
            let mut writer = Writer {
                writer: handle,
                linewrap: Some(4),
            };
            let record = Record::with_attrs("id", Some("desc"), b"ACGT");

            let write_result = writer.write_record(&record);
            assert!(write_result.is_ok());
        }

        let actual = fs::read_to_string(path).unwrap();
        let expected = ">id desc\nACGT\n";

        assert!(fs::remove_file(path).is_ok());
        assert_eq!(actual, expected)
    }

    #[test]
    fn test_write_with_linewrap() {
        let width = 4;
        let mut writer = Writer::new(Vec::new());
        writer.set_linewrap(Some(width));
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE_WIDTH);

        let mut writer = Writer::with_capacity(100, Vec::new());
        writer.set_linewrap(Some(width));
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE_WIDTH);

        let mut writer = Writer::from_bufwriter(std::io::BufWriter::with_capacity(100, Vec::new()));
        writer.set_linewrap(Some(width));
        writer.write("id", Some("desc"), b"ACCGTAGGCTGA").unwrap();
        writer.write("id2", None, b"ATTGTTGTTTTA").unwrap();
        writer.flush().unwrap();
        assert_eq!(writer.writer.get_ref(), &WRITE_FASTA_FILE_WIDTH);
    }
}

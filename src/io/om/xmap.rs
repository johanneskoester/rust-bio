//! Structs amd trait to read files in XMAP format.
//!
//! Note: This module only supports version v0.2 of the XMAP format for now.
//!
//! # Example
//!
//! ## Reader
//!
//! In this example, we parse a XMAP file, iterate over its records. and output
//! some statistics.
//!
//! ```rust
//! use bio::io::om::xmap;
//!
//! let path = "tests/resources/valid_input.xmap";
//! let mut reader = xmap::Reader::from_path(&path).unwrap();
//!
//! let mut num_aligns = 0;
//!
//! while let Some(Ok(record)) = reader.next() {
//!     num_aligns += 1;
//! }
//!
//! println!("Number of alignments: {}", num_aligns);
//! ```
//!
//! ## Container
//!
//! If feasible, we can also build a container over the content of a XMAP file.
//!
//! This can be advantageous if out-of-order access is preferred or if the
//! users only wants to fetch alignments spanning a specific region.
//!
//! ```rust
//! use bio::io::om::xmap;
//!
//! let path = "tests/resources/valid_input.xmap";
//! let mut container = xmap::Container::from_path(&path).unwrap();
//!
//! println!("XMAP header: {:?}", container.header());
//! ```

use anyhow::{bail, Context, Result};
use derive_new::new;
use getset::Getters;
use lazy_static::lazy_static;
use regex::Regex;
use std::collections::{BTreeMap, HashMap};
use std::fmt;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::rc::Rc;

use crate::io::om::common::Error;
use crate::io::om::common::NextToErr;
use crate::io::om::common::Orientation;
use crate::io::om::common::ParseToErr;

lazy_static! {
    /// Regex command to describe globally valid alignment strings in XMAP.
    static ref RE_ALIGN_GLOBAL: Regex =
        Regex::new(r"\A(\([0-9]+,[0-9]+\))*\z").unwrap();
    /// Regex command to describe locally valid alignment strings in XMAP.
    static ref RE_ALIGN_LOCAL: Regex =
        Regex::new(r"\((?<ref_id>[0-9]+),(?<qry_id>[0-9]+)\)").unwrap();
    /// Regex command to describe valid CIGAR strings in XMAP.
    static ref RE_CIGAR: Regex = Regex::new(r"\A([0-9]+[MID])*\z").unwrap();
}

/// A XMAP record.
#[derive(Clone, Debug, Getters, PartialEq)]
pub struct Record {
    #[getset(get = "pub")]
    id: u32,
    #[getset(get = "pub")]
    qry_id: u32,
    #[getset(get = "pub")]
    ref_id: u32,
    qry_start: f64,
    qry_end: f64,
    #[getset(get = "pub")]
    ref_start: f64,
    ref_end: f64,
    orientation: Orientation,
    confidence: u32,
    #[getset(get = "pub")]
    hit_enum: Vec<u8>,
    #[getset(get = "pub")]
    qry_len: f64,
    ref_len: f64,
    label_channel: u8,
    #[getset(get = "pub")]
    alignment: Vec<(u32, u32)>,
}

impl Record {
    /// Constructs a new `Record`.
    pub fn new() -> Record {
        Record {
            id: 0,
            qry_id: 0,
            ref_id: 0,
            qry_start: -1.0,
            qry_end: -1.0,
            ref_start: -1.0,
            ref_end: -1.0,
            orientation: Orientation::Unknown,
            confidence: 0,
            hit_enum: Vec::new(),
            qry_len: -1.0,
            ref_len: -1.0,
            label_channel: 0,
            alignment: Vec::new(),
        }
    }

    /// Fills based on given alignment line in a XMAP.
    ///
    /// # Example
    /// ```rust
    /// use bio::io::om::xmap::Record;
    ///
    /// let line = "1\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)";
    ///
    /// let mut rec = Record::new();
    /// rec.fill(&line);
    /// ```
    pub fn fill(&mut self, line: &str) -> Result<()> {
        let mut split = line.split('\t');

        self.id = split.try_parse("XmapEntryID")?;
        self.qry_id = split.try_parse("QryContigID")?;
        self.ref_id = split.try_parse("RefContigID")?;
        self.qry_start = split.try_parse("QryStartPos")?;
        self.qry_end = split.try_parse("QryEndPos")?;
        self.ref_start = split.try_parse("RefStartPos")?;
        self.ref_end = split.try_parse("RefEndPos")?;
        self.orientation = Record::read_orientation(split.try_next("Orientation")?)?;
        self.confidence = split.try_parse("Confidence")?;
        self.hit_enum = Record::is_valid_cigar(String::from(split.try_next("HitEnum")?))?;
        self.qry_len = split.try_parse("QryLen")?;
        self.ref_len = split.try_parse("RefLen")?;
        self.label_channel = split.try_parse("LabelChannel")?;
        self.alignment = Record::read_alignment(split.try_next("Alignment")?)?;

        Ok(())
    }

    /// Constructs a new `Record` based on alignment line of a XMAP.
    pub fn from(line: &str) -> Result<Self> {
        let mut rec = Record::new();
        let _ = rec.fill(&line)?;
        Ok(rec)
    }

    /// Reads orientation from string.
    fn read_orientation(dir: &str) -> Result<Orientation> {
        if dir == "+" {
            Ok(Orientation::Forward)
        } else if dir == "-" {
            Ok(Orientation::Backward)
        } else {
            bail!(Error::InvalidOrientation(String::from(dir)))
        }
    }

    /// Checks whether a CIGAR string is valid.
    fn is_valid_cigar(cigar: String) -> Result<Vec<u8>> {
        if RE_CIGAR.is_match(cigar.as_str()) {
            Ok(Vec::from(cigar))
        } else {
            bail!(Error::InvalidCigar)
        }
    }

    /// Reads alignment as vector of pairings between reference and query IDs.
    fn read_alignment(alignment: &str) -> Result<Vec<(u32, u32)>> {
        if RE_ALIGN_GLOBAL.is_match(alignment) {
            let align_vec: Vec<(u32, u32)> = RE_ALIGN_LOCAL
                .captures_iter(alignment)
                .map(|entry| {
                    (
                        entry["ref_id"].parse().unwrap(),
                        entry["qry_id"].parse().unwrap(),
                    )
                })
                .collect();
            Ok(align_vec)
        } else {
            bail!(Error::InvalidAlignment);
        }
    }
}

/// Trait for XMAP readers.
pub trait XmapRead: Iterator<Item = Result<Record>> {
    /// Reads data of next alignment in XMAP into the given [Record](struct.Record.html).
    /// An empty record indicates that no more records can be read.
    fn read_into(&mut self, record: &mut Record) -> Result<bool>;
}

/// A XMAP reader.
pub struct Reader<R: BufRead> {
    stream: R,
    header: Vec<String>,
    buffer: String,
}

impl Reader<BufReader<File>> {
    /// Reads XMAP from given file path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let stream = BufReader::new(File::open(path).context(Error::InvalidPath)?);
        Reader::from_stream(stream)
    }
}

impl<R: BufRead> Reader<R> {
    /// Reads XMAP from buffered stream.
    pub fn from_stream(mut stream: R) -> Result<Self> {
        let mut header = Vec::new();
        let mut buffer = String::new();
        loop {
            buffer.clear();
            if stream.read_line(&mut buffer)? == 0 {
                break;
            };
            if buffer.starts_with('#') {
                header.push(String::from(buffer.trim_end()));
            } else {
                break;
            }
        }
        Ok(Reader {
            stream,
            header,
            buffer,
        })
    }

    /// Constructs [Container](struct.Container.html) over XMAP content while consuming reader.
    ///
    /// # Example
    /// ```rust
    /// use bio::io::om::xmap::Reader;
    ///
    /// let path = "tests/resources/valid_input.xmap";
    /// let reader = Reader::from_path(&path).unwrap();
    /// let container = reader.into_container();
    /// ```
    pub fn into_container(self) -> Result<Container> {
        let header = self.header.clone();
        let mut inner = HashMap::new();
        for rec in self {
            let rec = match rec {
                Ok(rec) => rec,
                Err(e) => return Err(e),
            };
            inner
                .entry(*rec.ref_id())
                .or_insert_with(BTreeMap::new)
                .insert(*rec.ref_start() as u64, rec.into());
        }
        Ok(Container::new(header, inner))
    }
}

impl<R: BufRead> XmapRead for Reader<R> {
    fn read_into(&mut self, record: &mut Record) -> Result<bool> {
        if self.buffer.is_empty() {
            return Ok(false);
        }
        let res = match record.fill(self.buffer.trim()) {
            Ok(()) => Ok(true),
            Err(e) => Err(e),
        };
        self.buffer.clear();
        match self.stream.read_line(&mut self.buffer) {
            Ok(_) => res,
            Err(e) => Err(e.into()),
        }
    }
}

/// Iterator over records.
impl<R: BufRead> Iterator for Reader<R> {
    type Item = Result<Record>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut rec = Record::new();
        match self.read_into(&mut rec) {
            Ok(true) => Some(Ok(rec)),
            Ok(false) => None,
            Err(e) => Some(Err(e)),
        }
    }
}

/// A container with XMAP header and records stored in B-trees.
///
/// Each tree corresponds to one reference contig and is ordered by start position.
#[derive(Clone, Debug, Getters, new, PartialEq)]
pub struct Container {
    #[getset(get = "pub")]
    header: Vec<String>,
    pos_trees: HashMap<u32, BTreeMap<u64, Rc<Record>>>,
}

impl Container {
    /// Constructs a new `Container` from path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = Reader::from_path(path)?;
        reader.into_container()
    }

    /// Fetches an iterator over records contained in a given region in a contig.
    pub fn fetch(
        &self,
        region_contig: u32,
        region_start: u64,
        region_end: u64,
    ) -> Result<impl Iterator<Item = &Rc<Record>>> {
        let hits = match self.pos_trees.get(&region_contig) {
            Some(hit_tree) => hit_tree.range(region_start..region_end),
            None => Default::default(),
        };
        Ok(hits.map(|(_, rec)| rec))
    }
}

impl fmt::Display for Container {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        f.debug_struct("Container")
            .field("header", &format_args!("{:?}", &self.header))
            .field(
                "pos_trees",
                &format_args! {
                    "{:?}",
                    self.pos_trees.clone().iter_mut().map(
                        |(&cont, tree)|
                        (
                            cont, tree.iter_mut().map(
                                |(&id, rec)|
                                (id, rec.id())
                            ).collect::<HashMap<_, _>>()
                        )
                    ).collect::<HashMap<_, _>>()
                },
            )
            .finish()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_create_new_record() {
        assert_eq!(
            Record::new(),
            Record {
                id: 0,
                qry_id: 0,
                ref_id: 0,
                qry_start: -1.0,
                qry_end: -1.0,
                ref_start: -1.0,
                ref_end: -1.0,
                orientation: Orientation::Unknown,
                confidence: 0,
                hit_enum: Vec::new(),
                qry_len: -1.0,
                ref_len: -1.0,
                label_channel: 0,
                alignment: Vec::new()
            }
        );
    }

    #[test]
    fn test_read_forward_orientation() {
        assert_eq!(Record::read_orientation("+").unwrap(), Orientation::Forward,);
    }

    #[test]
    fn test_read_backward_orientation() {
        assert_eq!(
            Record::read_orientation("-").unwrap(),
            Orientation::Backward,
        );
    }

    #[test]
    fn test_read_invalid_orientation() {
        assert_eq!(
            format!("{}", Record::read_orientation("#").unwrap_err()),
            "InvalidData: Orientation # is not defined.",
        );
    }

    #[test]
    fn test_is_valid_cigar_with_valid_input() {
        assert_eq!(
            Record::is_valid_cigar(String::from("2M3D5I6M")).unwrap(),
            Vec::from("2M3D5I6M"),
        );
    }

    #[test]
    fn test_is_valid_cigar_without_count() {
        assert_eq!(
            format!("{}", Record::is_valid_cigar(String::from("M")).unwrap_err()),
            "InvalidData: HitEnum is an invalid CIGAR string.",
        );
    }

    #[test]
    fn test_is_valid_cigar_without_indicator() {
        assert_eq!(
            format!("{}", Record::is_valid_cigar(String::from("2")).unwrap_err()),
            "InvalidData: HitEnum is an invalid CIGAR string.",
        );
    }

    #[test]
    fn test_is_valid_cigar_with_invalid_indicator() {
        assert_eq!(
            format!(
                "{}",
                Record::is_valid_cigar(String::from("2MM")).unwrap_err()
            ),
            "InvalidData: HitEnum is an invalid CIGAR string.",
        );
    }

    #[test]
    fn test_read_alignment_with_valid_input() {
        let mut output_vec = Vec::new();
        output_vec.push((5, 6));
        output_vec.push((6, 7));
        output_vec.push((8, 8));

        assert_eq!(
            Record::read_alignment("(5,6)(6,7)(8,8)").unwrap(),
            output_vec,
        );
    }

    #[test]
    fn test_read_alignment_with_invalid_input() {
        assert_eq!(
            format!("{}", Record::read_alignment("(5,6)(,7)(8,f)6").unwrap_err()),
            "InvalidData: Alignment is not formatted correctly.",
        );
    }

    #[test]
    fn test_fill_record_with_valid_input() {
        let mut rec = Record::new();
        let line = String::from("1\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        assert!(rec.fill(&line).is_ok());

        let mut output_vec = Vec::new();
        output_vec.push((5, 3));
        output_vec.push((6, 2));

        assert_eq!(
            rec,
            Record {
                id: 1,
                qry_id: 2,
                ref_id: 3,
                qry_start: 4.1,
                qry_end: 5.3,
                ref_start: 6.0,
                ref_end: 7.0,
                orientation: Orientation::Backward,
                confidence: 8,
                hit_enum: Vec::from("2M"),
                qry_len: 1.2,
                ref_len: 1.0,
                label_channel: 1,
                alignment: output_vec
            }
        );
    }

    #[test]
    fn test_fill_record_with_invalid_id() {
        let mut rec = Record::new();
        let line = String::from("1f5\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: XmapEntryID is not a valid u32 integer."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_qry_id() {
        let mut rec = Record::new();
        let line = String::from(
            "15\t2.4\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\
            \t1\t(5,3)(6,2)",
        );
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: QryContigID is not a valid u32 integer."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_ref_id() {
        let mut rec = Record::new();
        let line = String::from(
            "15\t2\t3.2\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\
            \t1\t(5,3)(6,2)",
        );
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: RefContigID is not a valid u32 integer."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_qry_start() {
        let mut rec = Record::new();
        let line = String::from(
            "15\t2\t3\tf\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0
            \t1\t(5,3)(6,2)",
        );
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: QryStartPos is not a valid f64 float."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_qry_end() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\tf\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: QryEndPos is not a valid f64 float."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_ref_start() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\tf\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: RefStartPos is not a valid f64 float."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_ref_end() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\tf\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: RefEndPos is not a valid f64 float."
        );
    }

    #[test]
    fn test_fill_record_without_orientation() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Orientation."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_orientation() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t#\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: Orientation # is not defined."
        );
    }

    #[test]
    fn test_fill_record_without_confidence() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Confidence."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_confidence() {
        let mut rec = Record::new();
        let line = String::from(
            "15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8.5\t2M\t1.2\t1.0\
            \t1\t(5,3)(6,2)",
        );
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: Confidence is not a valid u32 integer."
        );
    }

    #[test]
    fn test_fill_record_without_hit_enum() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field HitEnum."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_hit_enum() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2\t1.2\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: HitEnum is an invalid CIGAR string."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_qry_len() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\tf\t1.0\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: QryLen is not a valid f64 float."
        );
    }

    #[test]
    fn test_fill_record_without_ref_len() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field RefLen."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_ref_len() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\tf\t1\t(5,3)(6,2)");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: RefLen is not a valid f64 float."
        );
    }

    #[test]
    fn test_fill_record_without_label_channel() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field LabelChannel."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_label_channel() {
        let mut rec = Record::new();
        let line = String::from(
            "15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\
            \t1.5\t(5,3)(6,2)",
        );
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: LabelChannel is not a valid u8 integer."
        );
    }

    #[test]
    fn test_fill_record_without_alignment() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Alignment."
        );
    }

    #[test]
    fn test_fill_record_with_invalid_alignment() {
        let mut rec = Record::new();
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)5");
        let result = rec.fill(&line);
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: Alignment is not formatted correctly."
        );
    }

    #[test]
    fn test_from_record_with_valid_input() {
        let line = String::from("1\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let rec = Record::from(&line);

        let mut output_vec = Vec::new();
        output_vec.push((5, 3));
        output_vec.push((6, 2));

        assert_eq!(
            rec.unwrap(),
            Record {
                id: 1,
                qry_id: 2,
                ref_id: 3,
                qry_start: 4.1,
                qry_end: 5.3,
                ref_start: 6.0,
                ref_end: 7.0,
                orientation: Orientation::Backward,
                confidence: 8,
                hit_enum: Vec::from("2M"),
                qry_len: 1.2,
                ref_len: 1.0,
                label_channel: 1,
                alignment: output_vec
            }
        );
    }

    #[test]
    fn test_from_record_with_invalid_input() {
        let line = String::from("1\t2\t3\t4.1\t5.3\t6.0\t7.0\t-");
        let rec = Record::from(&line);

        assert_eq!(
            format!("{}", rec.unwrap_err()),
            "Truncated file: Cannot extract field Confidence.",
        )
    }

    #[test]
    fn test_from_record_equals_fill_record() {
        let line = String::from(
            "36\t1\t17\t5119.30\t465.17\t3888.0\t9624.0\t-\t8606\
            \t4M1D2M1D1M1D1M1D1M\t4654.13\t5736.0\t1\
            \t(37,10)(38,9)(39,8)(40,7)(42,6)(43,5)(45,4)(47,3)(49,2)",
        );

        let mut rec_fill = Record::new();
        assert!(rec_fill.fill(&line).is_ok());

        let rec_from = Record::from(&line);
        assert!(rec_from.is_ok());

        assert_eq!(rec_fill, rec_from.unwrap());
    }

    #[test]
    fn test_get_id() {
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let rec = Record::from(&line).unwrap();

        assert_eq!(*rec.id(), 15)
    }

    #[test]
    fn test_get_ref_id() {
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let rec = Record::from(&line).unwrap();

        assert_eq!(*rec.ref_id(), 3)
    }

    #[test]
    fn test_get_ref_start() {
        let line = String::from("15\t2\t3\t4.1\t5.3\t6.0\t7.0\t-\t8\t2M\t1.2\t1.0\t1\t(5,3)(6,2)");
        let rec = Record::from(&line).unwrap();

        assert_eq!(*rec.ref_start(), 6.0)
    }

    #[test]
    fn test_from_path_container_with_valid_path() {
        let container = Container::from_path("tests/resources/valid_input.xmap");
        let mut header = Vec::new();
        header.push(String::from("# XMAP File Version:\t0.2"));
        header.push(String::from("# Label Channels:\t1"));
        header.push(String::from(
            "#h\tXmapEntryID\tQryContigID\tRefContigID\tQryStartPos\
             \tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\
             \tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment",
        ));
        header.push(String::from(
            "#f\tint\tint\tint\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\
             \tstring\tfloat\tfloat\tint\tstring",
        ));

        let mut tree_17 = BTreeMap::new();
        tree_17.insert(
            3888 as u64,
            Record::from(
                "36\t1\t17\t5119.30\t465.17\t3888.0\t9624.0\t-\t8606\
                \t4M1D2M1D1M1D1M1D1M\t4654.13\t5736.0\t1\
                \t(37,10)(38,9)(39,8)(40,7)(42,6)(43,5)(45,4)(47,3)(49,2)",
            )
            .unwrap()
            .into(),
        );
        tree_17.insert(
            5396 as u64,
            Record::from(
                "37\t2\t17\t2617.12\t5942.70\t5396.0\t8148.0\t+\t4434\
                \t1M1D1M1D6M\t3325.58\t2752.0\t1\
                \t(22,1)(24,2)(26,3)(27,4)(28,5)(29,6)(30,7)(31,8)",
            )
            .unwrap()
            .into(),
        );
        tree_17.insert(
            6152 as u64,
            Record::from(
                "38\t3\t17\t538.88\t4004.65\t6152.0\t1568.0\t+\t1726\
                 \t1M1D2M2D3M1D1M\t3465.77\t4594.0\t1\
                 \t(62,1)(64,2)(65,3)(68,4)(69,5)(70,6)(72,7)",
            )
            .unwrap()
            .into(),
        );
        let mut tree_22 = BTreeMap::new();
        tree_22.insert(
            4974 as u64,
            Record::from(
                "40\t4\t22\t477.62\t2933.80\t4974.0\t7257.0\t+\t3416\
                \t2M2D3M1D1M1D1M\t2456.18\t3283.0\t1\
                \t(84,1)(85,2)(88,3)(89,4)(90,5)(92,6)(94,7)",
            )
            .unwrap()
            .into(),
        );
        let mut tree_map = HashMap::new();
        tree_map.insert(22 as u32, tree_22);
        tree_map.insert(17 as u32, tree_17);

        assert_eq!(
            container.unwrap(),
            Container {
                header,
                pos_trees: tree_map
            },
        );
    }

    #[test]
    fn test_from_path_container_with_invalid_input() {
        let res = Container::from_path("tests/resources/invalid_input.xmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Truncated file: Cannot extract field RefContigID.",
        )
    }

    #[test]
    fn test_from_path_container_with_invalid_path() {
        let res = Container::from_path("/not/a/real/path.xmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Invalid Path: Cannot locate specified path.",
        )
    }

    #[test]
    fn test_from_path_container_with_empty_file() {
        let container = Container::from_path("tests/resources/empty_file.txt").unwrap();
        assert_eq!(
            container,
            Container {
                header: Vec::new(),
                pos_trees: HashMap::new(),
            },
        )
    }

    #[test]
    fn test_from_path_container_with_non_utf8_file_without_record() {
        let res = Container::from_path("tests/resources/non_utf8_without_record.txt");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "stream did not contain valid UTF-8",
        )
    }

    #[test]
    fn test_from_path_container_with_non_utf8_file_with_record() {
        let res = Container::from_path("tests/resources/non_utf8_with_record.xmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "stream did not contain valid UTF-8",
        )
    }

    #[test]
    fn test_fetch_with_existing_contig() {
        let container = Container::from_path("tests/resources/valid_input.xmap").unwrap();

        let mut trees = Vec::new();
        let rec_1 = Record::from(
            "36\t1\t17\t5119.30\t465.17\t3888.0\t9624.0\t-\t8606\
            \t4M1D2M1D1M1D1M1D1M\t4654.13\t5736.0\t1\
            \t(37,10)(38,9)(39,8)(40,7)(42,6)(43,5)(45,4)(47,3)(49,2)",
        )
        .unwrap()
        .into();
        let rec_2 = Record::from(
            "37\t2\t17\t2617.12\t5942.70\t5396.0\t8148.0\t+\t4434\
            \t1M1D1M1D6M\t3325.58\t2752.0\t1\
            \t(22,1)(24,2)(26,3)(27,4)(28,5)(29,6)(30,7)(31,8)",
        )
        .unwrap()
        .into();
        trees.push(&rec_1);
        trees.push(&rec_2);

        let mut range = container.fetch(17, 1, 70000000).unwrap();

        let mut count = 0;
        loop {
            assert_eq!(range.next().unwrap(), trees[count],);
            count += 1;
            if count == trees.len() {
                break;
            }
        }
    }

    #[test]
    fn test_fetch_without_existing_contig() {
        let container = Container::from_path("tests/resources/valid_input.xmap").unwrap();

        assert_eq!(container.fetch(1, 43, 234).unwrap().count(), 0)
    }

    #[test]
    fn test_fmt_display() {
        let container = Container::from_path("tests/resources/single_line.xmap").unwrap();
        assert_eq!(
            container.to_string(),
            "Container { header: [\"# XMAP File Version:\\t0.2\", \
            \"# Label Channels:\\t1\", \
            \"#h\\tXmapEntryID\\tQryContigID\\tRefContigID\\tQryStartPos\
            \\tQryEndPos\\tRefStartPos\\tRefEndPos\\tOrientation\\tConfidence\
            \\tHitEnum\\tQryLen\\tRefLen\\tLabelChannel\\tAlignment\", \
            \"#f\\tint\\tint\\tint\\tfloat\\tfloat\\tfloat\\tfloat\\tstring\
            \\tfloat\\tstring\\tfloat\\tfloat\\tint\\tstring\"], \
            pos_trees: {22: {4974: 40}} }"
        )
    }

    #[test]
    fn test_get_header() {
        let container = Container::from_path("tests/resources/valid_input.xmap").unwrap();
        let mut header = Vec::new();
        header.push(String::from("# XMAP File Version:\t0.2"));
        header.push(String::from("# Label Channels:\t1"));
        header.push(String::from(
            "#h\tXmapEntryID\tQryContigID\tRefContigID\tQryStartPos\
            \tQryEndPos\tRefStartPos\tRefEndPos\tOrientation\tConfidence\
            \tHitEnum\tQryLen\tRefLen\tLabelChannel\tAlignment",
        ));
        header.push(String::from(
            "#f\tint\tint\tint\tfloat\tfloat\tfloat\tfloat\tstring\tfloat\
            \tstring\tfloat\tfloat\tint\tstring",
        ));

        assert_eq!(container.header(), &header,)
    }
}

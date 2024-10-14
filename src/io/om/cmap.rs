//! Structs amd trait to read files in CMAP format.
//!
//! Note: This module only supports version v0.1 of the CMAP format for now.
//!
//! # Example
//!
//! ## Reader
//!
//! In this example, we parse a CMAP file, iterate over its records. and output
//! some statistics.
//!
//! ```rust
//! use bio::io::om::cmap;
//!
//! let path = "tests/resources/valid_input.cmap";
//! let mut reader = cmap::Reader::from_path(&path).unwrap();
//!
//! let mut num_contigs = 0;
//! let mut num_labels = 0;
//!
//! while let Some(Ok(mut record)) = reader.next() {
//!     num_contigs += 1;
//!     num_labels += record.labels().len() - 1;
//! }
//!
//! println!("Number of contigs: {}", num_contigs);
//! println!("Number of labels: {}", num_labels);
//! ```
//!
//! ## Container
//!
//! If feasible, we can also build a container over the content of a CMAP file.
//!
//! This can be advantageous if out-of-order access is preferred, which is
//! often the case for CMAPs.
//!
//! ```rust
//! use bio::io::om::cmap;
//!
//! let path = "tests/resources/valid_input.cmap";
//! let mut container = cmap::Container::from_path(&path).unwrap();
//!
//! println!("CMAP header: {:?}", container.header());
//! println!("CMAP label: {:?}", container.labels());
//! ```

use anyhow::{bail, Context, Result};
use derive_new::new;
use getset::Getters;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::rc::Rc;

use crate::io::om::common::Error;
use crate::io::om::common::NextToErr;
use crate::io::om::common::ParseToErr;

/// A label as given in a CMAP file.
#[derive(Clone, Debug, Getters, PartialEq)]
pub struct Label {
    id: u32,
    #[getset(get = "pub")]
    pos: f64,
    std_dev: f32,
    coverage: u32,
    occurrence: u32,
}

/// A CMAP record.
#[derive(Clone, Debug, Getters, PartialEq)]
pub struct Record {
    #[getset(get = "pub")]
    id: u32,
    len: f64,
    label_channel: u8,
    #[getset(get = "pub")]
    labels: Vec<Label>,
}

impl Record {
    /// Constructs a new `Record`.
    pub fn new() -> Record {
        Record {
            id: 0,
            len: -1.0,
            label_channel: 0,
            labels: Vec::new(),
        }
    }

    /// Fills based on all lines corresponding to the same contig in a CMAP.
    ///
    /// # Example
    /// ```rust
    /// use bio::io::om::cmap::Record;
    ///
    /// let mut lines = Vec::new();
    /// lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
    /// lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
    /// lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
    /// lines.push("1\t248936424.0\t3\t4\t0\t248936424.0\t1.0\t1\t0");
    ///
    /// let mut rec = Record::new();
    /// rec.fill(lines);
    /// ```
    pub fn fill(&mut self, lines: Vec<&str>) -> Result<()> {
        let split = &mut lines[0].split('\t');

        self.id = split.try_parse("CMapID")?;
        self.len = split.try_parse("ContigLength")?;
        split.try_next("NumSites")?;
        split.try_next("SiteID")?;
        self.label_channel = split.try_parse("LabelChannel")?;
        split.try_next("Position")?;
        split.try_next("StdDev")?;
        split.try_next("Coverage")?;
        split.try_next("Occurrence")?;

        self.labels = Record::read_labels(lines)?;

        Ok(())
    }

    /// Constructs a new `Record` based on all lines corresponding to the same contig in a CMAP.
    pub fn from(lines: Vec<&str>) -> Result<Self> {
        let mut rec = Record::new();
        let _ = rec.fill(lines)?;
        Ok(rec)
    }

    /// Reads label information for all channels.
    fn read_labels(label_lines: Vec<&str>) -> Result<Vec<Label>> {
        let mut labels = Vec::new();
        for idx in 0..label_lines.len() {
            let mut line_split = label_lines[idx].split('\t');
            let _: u32 = line_split.try_parse("CMapID")?;
            line_split.try_next("ContigLength")?;
            line_split.try_next("NumSites")?;
            let id = line_split.try_parse("SiteID")?;
            line_split.try_next("LabelChannel")?;
            let pos = line_split.try_parse("Position")?;
            let std_dev = line_split.try_parse("StdDev")?;
            let coverage = line_split.try_parse("Coverage")?;
            let occurrence = line_split.try_parse("Occurrence")?;
            labels.push(Label {
                id,
                pos,
                std_dev,
                coverage,
                occurrence,
            });
        }
        Ok(labels)
    }
}

/// Trait for CMAP readers.
pub trait CmapRead: Iterator<Item = Result<Record>> {
    /// Reads data of next contig in CMAP into the given [Record](struct.Record.html).
    /// An empty record indicates that no more records can be read.
    fn read_into(&mut self, record: &mut Record) -> Result<bool>;
}

/// A CMAP reader.
pub struct Reader<R: BufRead> {
    stream: R,
    header: Vec<String>,
    buffer: String,
}

impl Reader<BufReader<File>> {
    /// Reads CMAP from given file path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let stream = BufReader::new(File::open(path).context(Error::InvalidPath)?);
        Reader::from_stream(stream)
    }

    /// Constructs [Container](struct.Container.html) over CMAP content while consuming reader.
    ///
    /// # Example
    /// ```rust
    /// use bio::io::om::cmap::Reader;
    ///
    /// let path = "tests/resources/valid_input.cmap";
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
            inner.insert(*rec.id(), rec.into());
        }
        Ok(Container::new(header, inner))
    }
}

impl<R: BufRead> Reader<R> {
    /// Reads CMAP from buffered stream.
    pub fn from_stream(mut stream: R) -> Result<Self> {
        let mut header = Vec::new();
        let mut buffer = String::new();
        let mut num_channels = 0;
        loop {
            buffer.clear();
            if stream.read_line(&mut buffer)? == 0 {
                break;
            };
            if buffer.starts_with('#') {
                if buffer.starts_with("# Label Channels:\t1") {
                    num_channels = 1;
                }
                header.push(String::from(buffer.trim_end()));
            } else {
                break;
            }
        }
        loop {
            let res = stream.read_line(&mut buffer)?;
            if res == 0 {
                bail!(Error::IncompleteCmapRecord);
            }
            if buffer.ends_with("\t0")
                || buffer.ends_with("\t0\n")
                || buffer.ends_with("\t0\r\n")
                || buffer.is_empty()
            {
                break;
            }
        }
        if num_channels == 0 {
            bail!(Error::InvalidLabelChannel("1"))
        } else {
            Ok(Reader {
                stream,
                header,
                buffer,
            })
        }
    }
}

impl<R: BufRead> CmapRead for Reader<R> {
    fn read_into(&mut self, record: &mut Record) -> Result<bool> {
        if self.buffer.is_empty() {
            return Ok(false);
        }
        let buf_lines: Vec<&str> = self.buffer.trim().lines().collect();
        let mut res = match record.fill(buf_lines) {
            Ok(()) => Ok(true),
            Err(e) => Err(e),
        };
        self.buffer.clear();
        loop {
            res = match self.stream.read_line(&mut self.buffer) {
                Ok(0) => {
                    if self.buffer.is_empty() {
                        res
                    } else {
                        bail!(Error::IncompleteCmapRecord)
                    }
                }
                Ok(_) => res,
                Err(e) => Err(e.into()),
            };
            if self.buffer.ends_with("\t0")
                || self.buffer.ends_with("\t0\n")
                || self.buffer.ends_with("\t0\r\n")
                || self.buffer.is_empty()
            {
                break;
            }
        }
        res
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

/// A container with CMAP header and records.
#[derive(Debug, Getters, new, PartialEq)]
pub struct Container {
    #[getset(get = "pub")]
    header: Vec<String>,
    contigs: HashMap<u32, Rc<Record>>,
}

impl Container {
    /// Constructs a new `Container` from path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = Reader::from_path(path)?;
        reader.into_container()
    }

    /// Returns label positions (for each contig respectively).
    pub fn labels(&mut self) -> HashMap<u32, Vec<f64>> {
        self.contigs
            .clone()
            .iter_mut()
            .map(|(id, rec)| {
                (
                    *id,
                    rec.labels()
                        .iter()
                        .map(|label| label.pos)
                        .collect::<Vec<_>>(),
                )
            })
            .collect::<HashMap<_, _>>()
    }

    /// Returns record with given ID.
    pub fn record(&self, id: u32) -> Result<&Rc<Record>> {
        if !self.contigs.contains_key(&id) {
            bail!(Error::InvalidKeyAccess(id))
        }
        Ok(&self.contigs[&id])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_label_pos() {
        let label = Label {
            id: 5,
            pos: 12428.0,
            std_dev: 1.0,
            coverage: 1,
            occurrence: 1,
        };

        assert_eq!(*label.pos(), 12428.0)
    }

    #[test]
    fn test_create_new_record() {
        assert_eq!(
            Record::new(),
            Record {
                id: 0,
                len: -1.0,
                label_channel: 0,
                labels: Vec::new(),
            }
        );
    }

    #[test]
    fn test_read_labels_with_valid_input() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        assert_eq!(
            Record::read_labels(label_lines).unwrap(),
            Vec::from([
                Label {
                    id: 1,
                    pos: 4454.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 1
                },
                Label {
                    id: 2,
                    pos: 27579.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 1
                },
                Label {
                    id: 3,
                    pos: 98003.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 1
                },
                Label {
                    id: 4,
                    pos: 248936424.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 0
                },
            ]),
        );
    }

    #[test]
    fn test_read_labels_without_cmap_id() {
        let mut label_lines = Vec::new();
        label_lines.push("");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "InvalidData: CMapID is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_read_labels_without_contig_length() {
        let mut label_lines = Vec::new();
        label_lines.push("1");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field ContigLength.",
        );
    }

    #[test]
    fn test_read_labels_without_num_sites() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field NumSites.",
        );
    }

    #[test]
    fn test_read_labels_without_site_id() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field SiteID.",
        );
    }

    #[test]
    fn test_read_labels_without_label_channel() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field LabelChannel.",
        );
    }

    #[test]
    fn test_read_labels_without_position() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field Position.",
        );
    }

    #[test]
    fn test_read_labels_without_std_dev() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field StdDev.",
        );
    }

    #[test]
    fn test_read_labels_without_coverage() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field Coverage.",
        );
    }

    #[test]
    fn test_read_labels_without_occurrence() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "Truncated file: Cannot extract field Occurrence.",
        );
    }

    #[test]
    fn test_fill_record_with_valid_input() {
        let mut rec = Record::new();

        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        assert!(rec.fill(label_lines.clone()).is_ok());

        assert_eq!(
            rec,
            Record {
                id: 1,
                len: 248936424.0,
                label_channel: 1,
                labels: Vec::from([
                    Label {
                        id: 1,
                        pos: 4454.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1
                    },
                    Label {
                        id: 2,
                        pos: 27579.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1
                    },
                    Label {
                        id: 3,
                        pos: 98003.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1
                    },
                    Label {
                        id: 4,
                        pos: 248936424.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 0
                    },
                ]),
            },
        )
    }

    #[test]
    fn test_fill_record_without_cmap_id() {
        let mut label_lines = Vec::new();
        label_lines.push("");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: CMapID is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_fill_record_without_contig_length() {
        let mut label_lines = Vec::new();
        label_lines.push("1");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field ContigLength.",
        );
    }

    #[test]
    fn test_fill_record_without_num_sites() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field NumSites.",
        );
    }

    #[test]
    fn test_fill_record_without_site_id() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field SiteID.",
        );
    }

    #[test]
    fn test_fill_record_without_label_channel() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field LabelChannel.",
        );
    }

    #[test]
    fn test_fill_record_without_position() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Position.",
        );
    }

    #[test]
    fn test_fill_record_without_std_dev() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field StdDev.",
        );
    }

    #[test]
    fn test_fill_record_without_coverage() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Coverage.",
        );
    }

    #[test]
    fn test_fill_record_without_occurrence() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Occurrence.",
        );
    }

    #[test]
    fn test_fill_record_without_complete_labels() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1");

        let mut rec = Record::new();
        let result = rec.fill(label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field Occurrence.",
        );
    }

    #[test]
    fn test_from_record_with_valid_input() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        let rec = Record::from(label_lines);

        assert_eq!(
            rec.unwrap(),
            Record {
                id: 1,
                len: 248936424.0,
                label_channel: 1,
                labels: Vec::from([
                    Label {
                        id: 1,
                        pos: 4454.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1
                    },
                    Label {
                        id: 2,
                        pos: 27579.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1
                    },
                    Label {
                        id: 3,
                        pos: 98003.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1
                    },
                    Label {
                        id: 4,
                        pos: 248936424.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 0
                    },
                ]),
            },
        );
    }

    #[test]
    fn test_from_record_with_invalid_input() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1");

        let rec = Record::from(label_lines);

        assert_eq!(
            format!("{}", rec.unwrap_err()),
            "Truncated file: Cannot extract field Occurrence.",
        )
    }

    #[test]
    fn test_from_record_equals_fill_record() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        let mut rec_fill = Record::new();
        assert!(rec_fill.fill(label_lines.clone()).is_ok());

        let rec_from = Record::from(label_lines.clone());
        assert!(rec_from.is_ok());

        assert_eq!(rec_fill, rec_from.unwrap());
    }

    #[test]
    fn test_get_id() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        let rec = Record::from(label_lines).unwrap();
        assert_eq!(*rec.id(), 1)
    }

    #[test]
    fn test_get_labels() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        let rec = Record::from(label_lines).unwrap();
        assert_eq!(
            rec.labels(),
            &Vec::from([
                Label {
                    id: 1,
                    pos: 4454.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 1
                },
                Label {
                    id: 2,
                    pos: 27579.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 1
                },
                Label {
                    id: 3,
                    pos: 98003.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 1
                },
                Label {
                    id: 4,
                    pos: 248936424.0,
                    std_dev: 1.0,
                    coverage: 1,
                    occurrence: 0
                },
            ])
        )
    }

    #[test]
    fn test_from_path_container_with_valid_input() {
        let container = Container::from_path("tests/resources/valid_input.cmap");

        let mut header = Vec::new();
        header.push(String::from("# CMAP File Version:\t0.1"));
        header.push(String::from("# Label Channels:\t1"));
        header.push(String::from("# Nickase Recognition Site 1:\tCTTAAG"));
        header.push(String::from("# Number of Consensus Nanomaps:\t2"));
        header.push(String::from(
            "#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\
            \tPosition\tStdDev\tCoverage\tOccurrence",
        ));
        header.push(String::from(
            "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint",
        ));

        let mut contigs = HashMap::new();
        contigs.insert(
            1,
            Record {
                id: 1,
                len: 248936424.0,
                label_channel: 1,
                labels: Vec::from([
                    Label {
                        id: 1,
                        pos: 4454.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1,
                    },
                    Label {
                        id: 2,
                        pos: 27579.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1,
                    },
                    Label {
                        id: 3,
                        pos: 98003.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1,
                    },
                    Label {
                        id: 4,
                        pos: 248936424.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 0,
                    },
                ]),
            }
            .into(),
        );
        contigs.insert(
            2,
            Record {
                id: 2,
                len: 242173531.0,
                label_channel: 1,
                labels: Vec::from([
                    Label {
                        id: 1,
                        pos: 5925.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1,
                    },
                    Label {
                        id: 2,
                        pos: 12065.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 1,
                    },
                    Label {
                        id: 3,
                        pos: 242173531.0,
                        std_dev: 1.0,
                        coverage: 1,
                        occurrence: 0,
                    },
                ]),
            }
            .into(),
        );

        assert_eq!(container.unwrap(), Container { header, contigs },);
    }

    #[test]
    fn test_from_path_container_with_invalid_channel() {
        let res = Container::from_path("tests/resources/invalid_channel.cmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "InvalidFormat: No valid label channel number was given. \
            Allowed number of label channels: 1.",
        )
    }

    #[test]
    fn test_from_path_container_with_invalid_input() {
        let res = Container::from_path("tests/resources/invalid_input.cmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Truncated file: Cannot extract field Occurrence.",
        )
    }

    #[test]
    fn test_from_path_container_with_incomplete_input() {
        let res = Container::from_path("tests/resources/incomplete_input.cmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Truncated file: Cannot extract complete CMAP record. The last \
            line has to have a value of 0 for the field Occurrence.",
        )
    }

    #[test]
    fn test_from_path_container_with_invalid_path() {
        let res = Container::from_path("/not/a/real/path.cmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Invalid Path: Cannot locate specified path.",
        )
    }

    #[test]
    fn test_from_path_container_with_empty_file() {
        let res = Container::from_path("tests/resources/empty_file.txt");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Truncated file: Cannot extract complete CMAP record. The last \
            line has to have a value of 0 for the field Occurrence.",
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
    fn test_from_path_container_with_non_utf8_file_with_incomplete_record() {
        let res = Container::from_path("tests/resources/non_utf8_with_incomplete_record.cmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "stream did not contain valid UTF-8",
        )
    }

    #[test]
    fn test_from_path_container_with_non_utf8_file_with_complete_record() {
        let res = Container::from_path("tests/resources/non_utf8_with_complete_record.cmap");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "stream did not contain valid UTF-8",
        )
    }

    #[test]
    fn test_get_header() {
        let container = Container::from_path("tests/resources/valid_input.cmap").unwrap();

        let mut header = Vec::new();
        header.push(String::from("# CMAP File Version:\t0.1"));
        header.push(String::from("# Label Channels:\t1"));
        header.push(String::from("# Nickase Recognition Site 1:\tCTTAAG"));
        header.push(String::from("# Number of Consensus Nanomaps:\t2"));
        header.push(String::from(
            "#h CMapId\tContigLength\tNumSites\tSiteID\tLabelChannel\
            \tPosition\tStdDev\tCoverage\tOccurrence",
        ));
        header.push(String::from(
            "#f int\tfloat\tint\tint\tint\tfloat\tfloat\tint\tint",
        ));

        assert_eq!(container.header(), &header,)
    }

    #[test]
    fn test_get_record_with_valid_key() {
        let container = Container::from_path("tests/resources/valid_input.cmap").unwrap();
        let mut label_lines = Vec::new();
        label_lines.push("1\t248936424.0\t3\t1\t1\t4454.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t2\t1\t27579.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t3\t1\t98003.0\t1.0\t1\t1");
        label_lines.push("1\t248936424.0\t3\t4\t1\t248936424.0\t1.0\t1\t0");

        let rec = Record::from(label_lines).unwrap();
        assert_eq!(container.record(1).unwrap(), &rec.into())
    }

    #[test]
    fn test_get_record_with_invalid_key() {
        let container = Container::from_path("tests/resources/valid_input.cmap").unwrap();
        let res = container.record(17);
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "InvalidKeyAccess: Key 17 does not exist and cannot be accessed.",
        )
    }

    #[test]
    fn test_get_contigs() {
        let mut container = Container::from_path("tests/resources/valid_input.cmap").unwrap();

        let mut output = HashMap::new();
        output.insert(1, Vec::from([4454.0, 27579.0, 98003.0, 248936424.0]));
        output.insert(2, Vec::from([5925.0, 12065.0, 242173531.0]));

        assert_eq!(container.labels(), output,)
    }
}

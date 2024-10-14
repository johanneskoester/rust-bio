//! Structs amd trait to read files in BNX format.
//!
//! Note: This module only supports version v1.2 of the BNX format for now.
//!
//! # Example
//!
//! ## Reader
//!
//! In this example, we parse a BNX file, iterate over its records. and output
//! some statistics.
//!
//! ```rust
//! use bio::io::om::bnx;
//!
//! let path = "tests/resources/valid_input.bnx";
//! let mut reader = bnx::Reader::from_path(&path).unwrap();
//!
//! let mut num_molecules = 0;
//! let mut num_labels_channel_1 = 0;
//! let mut num_labels_channel_2 = 0;
//!
//! while let Some(Ok(mut record)) = reader.next() {
//!     num_molecules += 1;
//!     num_labels_channel_1 += record.channel_1().len() - 1;
//!     num_labels_channel_2 += record.channel_2().len() - 1;
//! }
//!
//! println!("Number of molecules: {}", num_molecules);
//! println!("Number of labels in channel 1: {}", num_labels_channel_1);
//! println!("Number of labels in channel 2: {}", num_labels_channel_2);
//! ```
//!
//! ## Container
//!
//! If feasible, we can also build a container over the content of a BNX file.
//!
//! This can be advantageous if the BNX is unsorted or if out-of-order access
//! is preferred (e.g. when analyzing alignments, which are often not in
//! order).
//!
//! ```rust
//! use bio::io::om::bnx;
//!
//! let path = "tests/resources/valid_input.bnx";
//! let mut container = bnx::Container::from_path(&path).unwrap();
//!
//! println!("BNX header: {:?}", container.header());
//! println!("BNX label: {:?}", container.labels());
//! ```

use anyhow::{bail, Context, Result};
use derive_new::new;
use getset::Getters;
use itertools::izip;
use lazy_static::lazy_static;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;
use std::rc::Rc;

use crate::io::om::common::Error;
use crate::io::om::common::NextToErr;
use crate::io::om::common::Orientation;
use crate::io::om::common::ParseToErr;

lazy_static! {
    /// Regex command to describe valid position lines in BNX.
    static ref RE_POS_LINE: Regex =
        Regex::new(r"\A[12]\t([0-9.\t])*\z").unwrap();
    /// Regex command to describe valid SNR lines in BNX.
    static ref RE_SNR_LINE: Regex =
        Regex::new(r"\AQX[12]1\t([0-9.\t])*\z").unwrap();
    /// Regex command to describe valid intensity lines in BNX.
    static ref RE_INTENSITY_LINE: Regex =
        Regex::new(r"\AQX[12]2\t([0-9.\t])*\z").unwrap();
}

/// A label as given in a BNX file.
#[derive(Clone, Debug, PartialEq)]
pub struct Label {
    pos: f64,
    signal_to_noise_ratio: f32,
    intensity: f32,
}

/// A BNX record.
#[derive(Clone, Debug, Getters, PartialEq)]
pub struct Record {
    #[getset(get = "pub")]
    id: u32,
    len: f64,
    avg_intensity: f32,
    signal_to_noise_ratio: f32,
    num_labels: u32,
    original_id: u32,
    scan_number: u8,
    scan_direction: Orientation,
    chip_id: Vec<u8>,
    flowcell: u8,
    run_id: u32,
    global_scan_number: u32,
    #[getset(get = "pub")]
    channel_1: Vec<Label>,
    #[getset(get = "pub")]
    channel_2: Vec<Label>,
}

impl Record {
    /// Constructs a new `Record`.
    pub fn new() -> Record {
        Record {
            id: 0,
            len: -1.0,
            avg_intensity: -1.0,
            signal_to_noise_ratio: -1.0,
            num_labels: 0,
            original_id: 0,
            scan_number: 0,
            scan_direction: Orientation::Unknown,
            chip_id: Vec::new(),
            flowcell: 0,
            run_id: 0,
            global_scan_number: 0,
            channel_1: Vec::new(),
            channel_2: Vec::new(),
        }
    }

    /// Fills based on meta line (indicated by `0`) and label lines of a BNX.
    ///
    /// # Example
    /// ```rust
    /// use bio::io::om::bnx::Record;
    ///
    /// let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t1";
    /// let mut label_lines = Vec::new();
    /// label_lines.push("1\t23.4\t42.6\t134.5");
    /// label_lines.push("QX11\t4.5\t7.3");
    /// label_lines.push("QX12\t0.4\t0.8");
    ///
    /// let mut rec = Record::new();
    /// rec.fill(&meta_line, label_lines);
    /// ```
    pub fn fill(&mut self, meta_line: &str, label_lines: Vec<&str>) -> Result<()> {
        let mut split = meta_line.split('\t');
        split.next();

        self.id = split.try_parse("MoleculeID")?;
        self.len = split.try_parse("Length")?;
        self.avg_intensity = split.try_parse("AvgIntensity")?;
        self.signal_to_noise_ratio = split.try_parse("SNR")?;
        self.num_labels = split.try_parse("NumberofLabels")?;
        self.original_id = split.try_parse("OriginalMoleculeId")?;
        self.scan_number = split.try_parse("ScanNumber")?;
        self.scan_direction = Record::read_orientation(split.try_next("ScanDirection")?)?;
        self.chip_id = Vec::from(split.try_next("ChipId")?);
        self.flowcell = split.try_parse("Flowcell")?;
        self.run_id = split.try_parse("RunId")?;
        self.global_scan_number = split.try_parse("GlobalScanNumber")?;
        (self.channel_1, self.channel_2) = Record::read_labels(label_lines)?;

        Ok(())
    }

    /// Constructs a new `Record` based on meta line (indicated by `0`) and label lines of a BNX.
    pub fn from(meta_line: &str, label_lines: Vec<&str>) -> Result<Self> {
        let mut rec = Record::new();
        let _ = rec.fill(&meta_line, label_lines)?;
        Ok(rec)
    }

    /// Reads orientation from string.
    fn read_orientation(dir: &str) -> Result<Orientation> {
        if dir == "0" {
            Ok(Orientation::Forward)
        } else if dir == "1" {
            Ok(Orientation::Backward)
        } else if dir == "-1" {
            Ok(Orientation::Unknown)
        } else {
            bail!(Error::InvalidOrientation(String::from(dir)))
        }
    }

    /// Reads label information for all channels.
    fn read_labels(label_lines: Vec<&str>) -> Result<(Vec<Label>, Vec<Label>)> {
        let num_lines = label_lines.len();
        if num_lines == 3 {
            let channel_1 =
                Record::read_label_for_channel(label_lines[0], label_lines[1], label_lines[2])?;
            Ok((channel_1, Vec::new()))
        } else if num_lines == 6 {
            let channel_1 =
                Record::read_label_for_channel(label_lines[0], label_lines[2], label_lines[4])?;
            let channel_2 =
                Record::read_label_for_channel(label_lines[1], label_lines[3], label_lines[5])?;
            Ok((channel_1, channel_2))
        } else {
            bail!(Error::InvalidBnxRecord(num_lines));
        }
    }

    /// Reads label information for one channel.
    fn read_label_for_channel(
        pos_line: &str,
        snr_line: &str,
        intensity_line: &str,
    ) -> Result<Vec<Label>> {
        let mut labels: Vec<Label> = Vec::new();
        let pos_line = match RE_POS_LINE.is_match(pos_line) {
            true => &pos_line[2..],
            false => {
                bail!(Error::InvalidLabelLine("LabelChannel"))
            }
        };
        let snr_line = match RE_SNR_LINE.is_match(snr_line) {
            true => &snr_line[5..],
            false => {
                bail!(Error::InvalidLabelLine("Label SNR"))
            }
        };
        let intensity_line = match RE_INTENSITY_LINE.is_match(intensity_line) {
            true => &intensity_line[5..],
            false => bail!(Error::InvalidLabelLine("Label Intensity")),
        };
        let pos_count = pos_line.matches('\t').count();
        let snr_count = snr_line.matches('\t').count();
        let intensity_count = intensity_line.matches('\t').count();
        if pos_count - 1 != snr_count || snr_count != intensity_count {
            bail!(Error::InvalidBnxLabel)
        }
        let snr_line = &(snr_line.to_string() + "\t0.0");
        let intensity_line = &(intensity_line.to_string() + "\t0.0");
        for entry in izip!(
            pos_line.split('\t'),
            snr_line.split('\t'),
            intensity_line.split('\t')
        ) {
            labels.push(Label {
                pos: entry
                    .0
                    .parse()
                    .context(Error::InvalidType("LabelPosition", "f64 float"))?,
                signal_to_noise_ratio: entry
                    .1
                    .parse()
                    .context(Error::InvalidType("QualityScore (SNR)", "f32 float"))?,
                intensity: entry
                    .2
                    .parse()
                    .context(Error::InvalidType("QualityScore (Intensity)", "f32 float"))?,
            });
        }
        Ok(labels)
    }
}

/// Trait for BNX readers.
pub trait BnxRead: Iterator<Item = Result<Record>> {
    /// Reads data of next molecule in BNX into the given [Record](struct.Record.html).
    /// An empty record indicates that no more records can be read.
    fn read_into(&mut self, record: &mut Record) -> Result<bool>;
}

/// A BNX reader.
pub struct Reader<R: BufRead> {
    stream: R,
    header: Vec<String>,
    num_channels: u8,
    buffer: String,
}

impl Reader<BufReader<File>> {
    /// Reads BNX from given file path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let stream = BufReader::new(File::open(path).context(Error::InvalidPath)?);
        Reader::from_stream(stream)
    }

    /// Constructs [Container](struct.Container.html) over BNX content while consuming reader.
    ///
    /// # Example
    /// ```rust
    /// use bio::io::om::bnx::Reader;
    ///
    /// let path = "tests/resources/valid_input.bnx";
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
    /// Reads BNX from buffered stream.
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
                } else if buffer.starts_with("# Label Channels:\t2") {
                    num_channels = 2;
                } else {
                };
                header.push(String::from(buffer.trim_end()));
            } else {
                break;
            }
        }
        for _ in 0..3 * num_channels {
            stream.read_line(&mut buffer)?;
        }
        if num_channels == 0 {
            bail!(Error::InvalidLabelChannel("1 or 2"))
        } else {
            Ok(Reader {
                stream,
                header,
                num_channels,
                buffer,
            })
        }
    }
}

impl<R: BufRead> BnxRead for Reader<R> {
    fn read_into(&mut self, record: &mut Record) -> Result<bool> {
        if self.buffer.is_empty() {
            return Ok(false);
        }
        let buf_lines: Vec<&str> = self.buffer.trim().lines().collect();
        let mut res = match record.fill(buf_lines[0], buf_lines[1..].to_vec()) {
            Ok(()) => Ok(true),
            Err(e) => Err(e),
        };
        self.buffer.clear();
        for _ in 0..3 * self.num_channels + 1 {
            res = match self.stream.read_line(&mut self.buffer) {
                Ok(_) => res,
                Err(e) => Err(e.into()),
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

/// A container with BNX header and records.
#[derive(Debug, new, Getters, PartialEq)]
pub struct Container {
    #[getset(get = "pub")]
    header: Vec<String>,
    molecules: HashMap<u32, Rc<Record>>,
}

impl Container {
    /// Constructs a new `Container` from path.
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let reader = Reader::from_path(path)?;
        reader.into_container()
    }

    /// Returns label positions (for each molecule respectively).
    pub fn labels(&mut self) -> HashMap<u32, (Vec<f64>, Vec<f64>)> {
        self.molecules
            .clone()
            .iter_mut()
            .map(|(id, rec)| {
                (
                    *id,
                    (
                        rec.channel_1()
                            .iter()
                            .map(|label| label.pos)
                            .collect::<Vec<_>>(),
                        rec.channel_2()
                            .iter()
                            .map(|label| label.pos)
                            .collect::<Vec<_>>(),
                    ),
                )
            })
            .collect::<HashMap<_, _>>()
    }

    /// Returns record with given ID.
    pub fn record(&self, id: u32) -> Result<&Rc<Record>> {
        if !self.molecules.contains_key(&id) {
            bail!(Error::InvalidKeyAccess(id))
        }
        Ok(&self.molecules[&id])
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
                len: -1.0,
                avg_intensity: -1.0,
                signal_to_noise_ratio: -1.0,
                num_labels: 0,
                original_id: 0,
                scan_number: 0,
                scan_direction: Orientation::Unknown,
                chip_id: Vec::new(),
                flowcell: 0,
                run_id: 0,
                global_scan_number: 0,
                channel_1: Vec::new(),
                channel_2: Vec::new(),
            }
        );
    }

    #[test]
    fn test_read_forward_orientation() {
        assert_eq!(Record::read_orientation("0").unwrap(), Orientation::Forward,);
    }

    #[test]
    fn test_read_backward_orientation() {
        assert_eq!(
            Record::read_orientation("1").unwrap(),
            Orientation::Backward,
        );
    }

    #[test]
    fn test_read_unknown_orientation() {
        assert_eq!(
            Record::read_orientation("-1").unwrap(),
            Orientation::Unknown,
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
    fn test_read_label_for_channel_with_valid_input() {
        let pos_line = "1\t23.4\t42.6\t134.5";
        let snr_line = "QX11\t0.4\t0.8";
        let intensity_line = "QX12\t4.5\t7.3";

        let mut output = Vec::new();
        output.push(Label {
            pos: 23.4,
            signal_to_noise_ratio: 0.4,
            intensity: 4.5,
        });
        output.push(Label {
            pos: 42.6,
            signal_to_noise_ratio: 0.8,
            intensity: 7.3,
        });
        output.push(Label {
            pos: 134.5,
            signal_to_noise_ratio: 0.0,
            intensity: 0.0,
        });

        assert_eq!(
            Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap(),
            output,
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_label_format() {
        let pos_line = "1\t23.4\t42.6\t134.5";
        let snr_line = "QX11\t0.4\t0.8\t0.3";
        let intensity_line = "QX12\t4.5\t7.3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidFormat: Number of entries in label lines is not \
            compatible with BNX structure (one position more than both SNR \
            and intensity)."
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_label_channel_format() {
        let pos_line = "1\t23.4\t42.6\t13f4.5";
        let snr_line = "QX11\t0.4\t0.8";
        let intensity_line = "QX12\t4.5\t7.3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidData: LabelChannel is not formatted correctly.",
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_snr_channel_format() {
        let pos_line = "1\t23.4\t42.6\t134.5";
        let snr_line = "QX11\t0.4\tf";
        let intensity_line = "QX12\t4.5\t7.3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidData: Label SNR is not formatted correctly.",
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_intensity_channel_format() {
        let pos_line = "1\t23.4\t42.6\t134.5";
        let snr_line = "QX11\t0.4\t0.8";
        let intensity_line = "QX12\t4.5\t7f.3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidData: Label Intensity is not formatted correctly.",
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_label_channel_content() {
        let pos_line = "1\t23.4\t42.6\t..";
        let snr_line = "QX11\t0.4\t0.8";
        let intensity_line = "QX12\t4.5\t7.3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidData: LabelPosition is not a valid f64 float.",
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_snr_channel_content() {
        let pos_line = "1\t23.4\t42.6\t134.5";
        let snr_line = "QX11\t0.4\t.";
        let intensity_line = "QX12\t4.5\t7.3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidData: QualityScore (SNR) is not a valid f32 float.",
        )
    }

    #[test]
    fn test_read_label_for_channel_with_incorrect_intensity_channel_content() {
        let pos_line = "1\t23.4\t42.6\t134.5";
        let snr_line = "QX11\t0.4\t0.8";
        let intensity_line = "QX12\t4.5\t7..3";

        assert_eq!(
            format!(
                "{}",
                Record::read_label_for_channel(pos_line, snr_line, intensity_line).unwrap_err()
            ),
            "InvalidData: QualityScore (Intensity) is not a valid f32 float.",
        )
    }

    #[test]
    fn test_read_labels_with_valid_single_channel_data() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let mut output = Vec::new();
        output.push(Label {
            pos: 23.4,
            signal_to_noise_ratio: 0.4,
            intensity: 4.5,
        });
        output.push(Label {
            pos: 42.6,
            signal_to_noise_ratio: 0.8,
            intensity: 7.3,
        });
        output.push(Label {
            pos: 134.5,
            signal_to_noise_ratio: 0.0,
            intensity: 0.0,
        });

        assert_eq!(
            Record::read_labels(label_lines).unwrap(),
            (output, Vec::new()),
        )
    }

    #[test]
    fn test_read_labels_with_valid_two_channel_data() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("2\t23.5\t42.7\t134.6");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX21\t0.5\t0.9");
        label_lines.push("QX12\t4.5\t7.3");
        label_lines.push("QX22\t4.6\t7.4");

        let mut output_1 = Vec::new();
        output_1.push(Label {
            pos: 23.4,
            signal_to_noise_ratio: 0.4,
            intensity: 4.5,
        });
        output_1.push(Label {
            pos: 42.6,
            signal_to_noise_ratio: 0.8,
            intensity: 7.3,
        });
        output_1.push(Label {
            pos: 134.5,
            signal_to_noise_ratio: 0.0,
            intensity: 0.0,
        });

        let mut output_2 = Vec::new();
        output_2.push(Label {
            pos: 23.5,
            signal_to_noise_ratio: 0.5,
            intensity: 4.6,
        });
        output_2.push(Label {
            pos: 42.7,
            signal_to_noise_ratio: 0.9,
            intensity: 7.4,
        });
        output_2.push(Label {
            pos: 134.6,
            signal_to_noise_ratio: 0.0,
            intensity: 0.0,
        });

        assert_eq!(
            Record::read_labels(label_lines).unwrap(),
            (output_1, output_2),
        )
    }

    #[test]
    fn test_read_labels_with_invalid_single_channel_data() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t13f4.5");
        label_lines.push("QX11\t0.4");
        label_lines.push("QX12\t4..5");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "InvalidData: LabelChannel is not formatted correctly.",
        )
    }

    #[test]
    fn test_read_labels_with_invalid_two_channel_data_in_channel_1() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t13f4.5");
        label_lines.push("2\t23.5\t42.7\t134.6");
        label_lines.push("QX11\t0.4");
        label_lines.push("QX21\t0.5");
        label_lines.push("QX12\t4.5");
        label_lines.push("QX22\t4.6");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "InvalidData: LabelChannel is not formatted correctly.",
        )
    }

    #[test]
    fn test_read_labels_with_invalid_two_channel_data_in_channel_2() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("2\t23.5\t42.7\t13f4.6");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX21\t0.5\t0.9");
        label_lines.push("QX12\t4.5\t7.3");
        label_lines.push("QX22\t4.6\t7.4");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "InvalidData: LabelChannel is not formatted correctly.",
        )
    }

    #[test]
    fn test_read_labels_with_invalid_line_amount() {
        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t13f4.5");
        label_lines.push("2\t23.5\t42.7\t134.6");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX21\t0.5\t0.9");
        label_lines.push("QX12\t4.5\t7.3");

        assert_eq!(
            format!("{}", Record::read_labels(label_lines).unwrap_err()),
            "InvalidFormat: Number of lines (i.e. 5) is not compatible with \
            BNX record structure (only 4 or 7 are allowed).",
        )
    }

    #[test]
    fn test_fill_record_with_valid_input() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        assert!(rec.fill(&meta_line, label_lines.clone()).is_ok());

        let mut output = Vec::new();
        output.push(Label {
            pos: 23.4,
            signal_to_noise_ratio: 0.4,
            intensity: 4.5,
        });
        output.push(Label {
            pos: 42.6,
            signal_to_noise_ratio: 0.8,
            intensity: 7.3,
        });
        output.push(Label {
            pos: 134.5,
            signal_to_noise_ratio: 0.0,
            intensity: 0.0,
        });

        assert_eq!(
            rec,
            Record {
                id: 1,
                len: 134.5,
                avg_intensity: 0.6,
                signal_to_noise_ratio: 5.3,
                num_labels: 2,
                original_id: 1,
                scan_number: 3,
                scan_direction: Orientation::Unknown,
                chip_id: Vec::from("TBA"),
                flowcell: 4,
                run_id: 5,
                global_scan_number: 6,
                channel_1: output,
                channel_2: Vec::new(),
            }
        );
    }

    #[test]
    fn test_fill_record_with_invalid_id() {
        let mut rec = Record::new();
        let meta_line = "0\t1.8\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: MoleculeID is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_len() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t13f4.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: Length is not a valid f64 float.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_avg_intensity() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0f.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: AvgIntensity is not a valid f32 float.",
        );
    }

    #[test]
    fn test_fill_record_without_snr() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field SNR.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_snr() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5f.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: SNR is not a valid f32 float.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_num_labels() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2.1\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: NumberofLabels is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_original_id() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1.1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: OriginalMoleculeId is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_scan_number() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3.1\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: ScanNumber is not a valid u8 integer.",
        );
    }

    #[test]
    fn test_fill_record_without_scan_direction() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field ScanDirection.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_scan_direction() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t+\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: Orientation + is not defined.",
        );
    }

    #[test]
    fn test_fill_record_without_chip_id() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "Truncated file: Cannot extract field ChipId.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_flowcell() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4.1\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: Flowcell is not a valid u8 integer.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_run_id() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5.1\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: RunId is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_global_scan_number() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6.1";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: GlobalScanNumber is not a valid u32 integer.",
        );
    }

    #[test]
    fn test_fill_record_with_invalid_labels() {
        let mut rec = Record::new();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42f.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let result = rec.fill(&meta_line, label_lines.clone());
        assert!(result.is_err());

        assert_eq!(
            format!("{}", result.unwrap_err()),
            "InvalidData: LabelChannel is not formatted correctly.",
        );
    }

    #[test]
    fn test_from_record_with_valid_input() {
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let rec = Record::from(&meta_line, label_lines);

        let mut output = Vec::new();
        output.push(Label {
            pos: 23.4,
            signal_to_noise_ratio: 0.4,
            intensity: 4.5,
        });
        output.push(Label {
            pos: 42.6,
            signal_to_noise_ratio: 0.8,
            intensity: 7.3,
        });
        output.push(Label {
            pos: 134.5,
            signal_to_noise_ratio: 0.0,
            intensity: 0.0,
        });

        assert_eq!(
            rec.unwrap(),
            Record {
                id: 1,
                len: 134.5,
                avg_intensity: 0.6,
                signal_to_noise_ratio: 5.3,
                num_labels: 2,
                original_id: 1,
                scan_number: 3,
                scan_direction: Orientation::Unknown,
                chip_id: Vec::from("TBA"),
                flowcell: 4,
                run_id: 5,
                global_scan_number: 6,
                channel_1: output,
                channel_2: Vec::new(),
            }
        );
    }

    #[test]
    fn test_from_record_with_invalid_input() {
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let rec = Record::from(&meta_line, label_lines);

        assert_eq!(
            format!("{}", rec.unwrap_err()),
            "Truncated file: Cannot extract field Flowcell.",
        )
    }

    #[test]
    fn test_from_record_equals_fill_record() {
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let mut rec_fill = Record::new();
        assert!(rec_fill.fill(&meta_line, label_lines.clone()).is_ok());

        let rec_from = Record::from(&meta_line, label_lines.clone());
        assert!(rec_from.is_ok());

        assert_eq!(rec_fill, rec_from.unwrap());
    }

    #[test]
    fn test_get_id() {
        let meta_line = "0\t15\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let rec = Record::from(&meta_line, label_lines).unwrap();
        assert_eq!(*rec.id(), 15)
    }

    #[test]
    fn test_get_channel_1() {
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX12\t4.5\t7.3");

        let rec = Record::from(&meta_line, label_lines).unwrap();
        assert_eq!(
            rec.channel_1(),
            &Vec::from([
                Label {
                    pos: 23.4,
                    signal_to_noise_ratio: 0.4,
                    intensity: 4.5
                },
                Label {
                    pos: 42.6,
                    signal_to_noise_ratio: 0.8,
                    intensity: 7.3
                },
                Label {
                    pos: 134.5,
                    signal_to_noise_ratio: 0.0,
                    intensity: 0.0
                },
            ])
        )
    }

    #[test]
    fn test_get_channel_2() {
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t6";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("2\t23.5\t42.7\t134.6");
        label_lines.push("QX11\t0.4\t0.8");
        label_lines.push("QX21\t0.5\t0.9");
        label_lines.push("QX12\t4.5\t7.3");
        label_lines.push("QX22\t4.6\t7.4");

        let rec = Record::from(&meta_line, label_lines).unwrap();
        assert_eq!(
            rec.channel_2(),
            &Vec::from([
                Label {
                    pos: 23.5,
                    signal_to_noise_ratio: 0.5,
                    intensity: 4.6
                },
                Label {
                    pos: 42.7,
                    signal_to_noise_ratio: 0.9,
                    intensity: 7.4
                },
                Label {
                    pos: 134.6,
                    signal_to_noise_ratio: 0.0,
                    intensity: 0.0
                },
            ])
        )
    }

    #[test]
    fn test_from_path_container_with_valid_input() {
        println!("{:?}", std::env::current_dir());
        let container = Container::from_path("tests/resources/valid_input.bnx");

        let mut header = Vec::new();
        header.push(String::from("# BNX File Version:\t1.2"));
        header.push(String::from("# Label Channels:\t2"));
        header.push(String::from("# Nickase Recognition Site 1:\tCTTAAG"));

        let mut molecules = HashMap::new();
        molecules.insert(
            1,
            Record {
                id: 1,
                len: 134.5,
                avg_intensity: 0.6,
                signal_to_noise_ratio: 5.3,
                num_labels: 2,
                original_id: 1,
                scan_number: 3,
                scan_direction: Orientation::Unknown,
                chip_id: Vec::from("TBA"),
                flowcell: 4,
                run_id: 5,
                global_scan_number: 1,
                channel_1: Vec::from([
                    Label {
                        pos: 23.4,
                        signal_to_noise_ratio: 4.5,
                        intensity: 0.4,
                    },
                    Label {
                        pos: 42.6,
                        signal_to_noise_ratio: 7.3,
                        intensity: 0.8,
                    },
                    Label {
                        pos: 134.5,
                        signal_to_noise_ratio: 0.0,
                        intensity: 0.0,
                    },
                ]),
                channel_2: Vec::from([
                    Label {
                        pos: 23.5,
                        signal_to_noise_ratio: 4.6,
                        intensity: 0.5,
                    },
                    Label {
                        pos: 42.7,
                        signal_to_noise_ratio: 7.4,
                        intensity: 0.9,
                    },
                    Label {
                        pos: 134.4,
                        signal_to_noise_ratio: 0.0,
                        intensity: 0.0,
                    },
                ]),
            }
            .into(),
        );
        molecules.insert(
            2,
            Record {
                id: 2,
                len: 350.87,
                avg_intensity: 0.15,
                signal_to_noise_ratio: 11.22,
                num_labels: 4,
                original_id: 2,
                scan_number: 1,
                scan_direction: Orientation::Unknown,
                chip_id: Vec::from("TBA"),
                flowcell: 1,
                run_id: 1,
                global_scan_number: 1,
                channel_1: Vec::from([
                    Label {
                        pos: 98.92,
                        signal_to_noise_ratio: 10.0380,
                        intensity: 0.0378,
                    },
                    Label {
                        pos: 120.27,
                        signal_to_noise_ratio: 8.1317,
                        intensity: 0.1308,
                    },
                    Label {
                        pos: 241.73,
                        signal_to_noise_ratio: 21.5573,
                        intensity: 0.0716,
                    },
                    Label {
                        pos: 252.47,
                        signal_to_noise_ratio: 9.8626,
                        intensity: 0.1518,
                    },
                    Label {
                        pos: 350.87,
                        signal_to_noise_ratio: 0.0,
                        intensity: 0.0,
                    },
                ]),
                channel_2: Vec::from([
                    Label {
                        pos: 98.8,
                        signal_to_noise_ratio: 9.0380,
                        intensity: 1.0378,
                    },
                    Label {
                        pos: 120.1,
                        signal_to_noise_ratio: 9.1317,
                        intensity: 1.1308,
                    },
                    Label {
                        pos: 241.6,
                        signal_to_noise_ratio: 20.5573,
                        intensity: 1.0716,
                    },
                    Label {
                        pos: 252.3,
                        signal_to_noise_ratio: 10.8626,
                        intensity: 1.1518,
                    },
                    Label {
                        pos: 350.7,
                        signal_to_noise_ratio: 0.0,
                        intensity: 0.0,
                    },
                ]),
            }
            .into(),
        );

        assert_eq!(container.unwrap(), Container { header, molecules },);
    }

    #[test]
    fn test_from_path_container_with_invalid_input() {
        let res = Container::from_path("tests/resources/invalid_input.bnx");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "Truncated file: Cannot extract field MoleculeID.",
        )
    }

    #[test]
    fn test_from_path_container_with_invalid_path() {
        let res = Container::from_path("/not/a/real/path.bnx");
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
            "InvalidFormat: No valid label channel number was given. \
            Allowed number of label channels: 1 or 2.",
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
        let res = Container::from_path("tests/resources/non_utf8_with_incomplete_record.bnx");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "stream did not contain valid UTF-8",
        )
    }

    #[test]
    fn test_from_path_container_with_non_utf8_file_with_complete_record() {
        let res = Container::from_path("tests/resources/non_utf8_with_complete_record.bnx");
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "stream did not contain valid UTF-8",
        )
    }

    #[test]
    fn test_get_header() {
        let container = Container::from_path("tests/resources/valid_input.bnx").unwrap();

        let mut header = Vec::new();
        header.push(String::from("# BNX File Version:\t1.2"));
        header.push(String::from("# Label Channels:\t2"));
        header.push(String::from("# Nickase Recognition Site 1:\tCTTAAG"));

        assert_eq!(container.header(), &header,)
    }

    #[test]
    fn test_get_record_with_valid_key() {
        let container = Container::from_path("tests/resources/valid_input.bnx").unwrap();
        let meta_line = "0\t1\t134.5\t0.6\t5.3\t2\t1\t3\t-1\tTBA\t4\t5\t1";

        let mut label_lines = Vec::new();
        label_lines.push("1\t23.4\t42.6\t134.5");
        label_lines.push("2\t23.5\t42.7\t134.4");
        label_lines.push("QX11\t4.5\t7.3");
        label_lines.push("QX21\t4.6\t7.4");
        label_lines.push("QX12\t0.4\t0.8");
        label_lines.push("QX22\t0.5\t0.9");

        let rec = Record::from(&meta_line, label_lines).unwrap();
        assert_eq!(container.record(1).unwrap(), &rec.into())
    }

    #[test]
    fn test_get_record_with_invalid_key() {
        let container = Container::from_path("tests/resources/valid_input.bnx").unwrap();
        let res = container.record(17);
        assert_eq!(
            format!("{}", res.unwrap_err()),
            "InvalidKeyAccess: Key 17 does not exist and cannot be accessed.",
        )
    }

    #[test]
    fn test_get_labels() {
        let mut container = Container::from_path("tests/resources/valid_input.bnx").unwrap();

        let mut output = HashMap::new();
        output.insert(
            1,
            (
                Vec::from([23.4, 42.6, 134.5]),
                Vec::from([23.5, 42.7, 134.4]),
            ),
        );
        output.insert(
            2,
            (
                Vec::from([98.92, 120.27, 241.73, 252.47, 350.87]),
                Vec::from([98.8, 120.1, 241.6, 252.3, 350.7]),
            ),
        );

        assert_eq!(container.labels(), output,)
    }
}

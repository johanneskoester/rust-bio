//! Common utilities for optical mapping data.

use anyhow::{Context, Result};

/// Enum to define possible alignment directions.
///
/// # Variants
/// * `Forward` - Forward direction (represented by `+` in XMAP and `0` in BNX).
/// * `Backward` - Backward direction (represented by `-` in XMAP and `1` in BNX).
/// * `Unknown` - Unknown direction (represented by `-1` in BNX).
#[derive(Clone, Debug, PartialEq)]
pub enum Orientation {
    Forward,
    Backward,
    Unknown,
}

/// Enum to define custom errors for optical mapping files.
///
/// # Variants
/// * `InvalidPath` - Path cannot be located.
/// * `IncompleteRecord` - Missing field in record.
/// * `InvalidLabelChannel` - No supported label channel was given.
/// * `InvalidBnxRecord` - Too few or many lines per BNX record.
/// * `InvalidBnxLabel` - Label lines in BNX include incomplete labels.
/// * `IncompleteCmapRecord` - Last record in CMAP is incomplete (Occurrence != 0).
/// * `InvalidLabelLine` - At least one label line in BNX is formatted incorrectly.
/// * `InvalidType` - Type of a field does not match expected value.
/// * `InvalidOrientation` - Field `ScanDirection` (for BNX) or `Orientation` (for CMAP) cannot be converted into valid [Orientation](enum.Orientation.html).
/// * `InvalidCigar` - Field `HitEnum` in XMAP does not contain a valid CIGAR.
/// * `InvalidAlignment` - Field `Alignment` in XMAP is formatted incorrectly.
/// * `InvalidKeyAccess` - Key is invalid and cannot be accessed.
#[derive(Clone, Debug, Eq, thiserror::Error, PartialEq)]
pub enum Error<'a> {
    #[error("Invalid Path: Cannot locate specified path.")]
    InvalidPath,

    #[error("Truncated file: Cannot extract field {0}.")]
    IncompleteRecord(String),

    #[error(
        "InvalidFormat: No valid label channel number was given. Allowed \
        number of label channels: {0}."
    )]
    InvalidLabelChannel(&'a str),

    #[error(
        "InvalidFormat: Number of lines (i.e. {0}) is not compatible with BNX \
        record structure (only 4 or 7 are allowed)."
    )]
    InvalidBnxRecord(usize),

    #[error(
        "InvalidFormat: Number of entries in label lines is not compatible \
        with BNX structure (one position more than both SNR and intensity)."
    )]
    InvalidBnxLabel,

    #[error(
        "Truncated file: Cannot extract complete CMAP record. The last line \
        has to have a value of 0 for the field Occurrence."
    )]
    IncompleteCmapRecord,

    #[error("InvalidData: {0} is not formatted correctly.")]
    InvalidLabelLine(&'a str),

    #[error("InvalidData: {0} is not a valid {1}.")]
    InvalidType(&'a str, &'a str),

    #[error("InvalidData: Orientation {0} is not defined.")]
    InvalidOrientation(String),

    #[error("InvalidData: HitEnum is an invalid CIGAR string.")]
    InvalidCigar,

    #[error("InvalidData: Alignment is not formatted correctly.")]
    InvalidAlignment,

    #[error("InvalidKeyAccess: Key {0} does not exist and cannot be accessed.")]
    InvalidKeyAccess(u32),
}

/// A helper trait to get next element from a split string.
pub trait NextToErr<'a> {
    fn try_next(&mut self, field: &'static str) -> Result<&'a str>;
}

impl<'a> NextToErr<'a> for std::str::Split<'a, char> {
    fn try_next(&mut self, field: &'static str) -> Result<&'a str> {
        self.next()
            .ok_or(Error::IncompleteRecord(String::from(field)).into())
    }
}

/// A helper trait to parse fields with numerical datatypes.
pub trait ParseToErr<'a, T> {
    fn try_parse(&mut self, field: &'static str) -> Result<T>;
}

impl<'a> ParseToErr<'a, u8> for std::str::Split<'a, char> {
    fn try_parse(&mut self, field: &'static str) -> Result<u8> {
        self.try_next(field)?
            .parse()
            .context(Error::InvalidType(field, "u8 integer"))
    }
}

impl<'a> ParseToErr<'a, u32> for std::str::Split<'a, char> {
    fn try_parse(&mut self, field: &'static str) -> Result<u32> {
        self.try_next(field)?
            .parse()
            .context(Error::InvalidType(field, "u32 integer"))
    }
}

impl<'a> ParseToErr<'a, f32> for std::str::Split<'a, char> {
    fn try_parse(&mut self, field: &'static str) -> Result<f32> {
        self.try_next(field)?
            .parse()
            .context(Error::InvalidType(field, "f32 float"))
    }
}

impl<'a> ParseToErr<'a, f64> for std::str::Split<'a, char> {
    fn try_parse(&mut self, field: &'static str) -> Result<f64> {
        self.try_next(field)?
            .parse()
            .context(Error::InvalidType(field, "f64 float"))
    }
}

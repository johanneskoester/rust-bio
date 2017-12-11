// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Various alignment and distance computing algorithms.

use utils::TextSlice;

pub mod pairwise;
pub mod distance;
pub mod sparse;


/// Alignment operations supported are match, substitution, insertion, deletion
/// and clipping. Clipping is a special boundary condition where you are allowed
/// to clip off the beginning/end of the sequence for a fixed clip penalty. The
/// clip penalty could be different for the two sequences x and y, and the
/// clipping operations on both are distinguishable (Xclip and Yclip). The usize
/// value associated with the clipping operations are the lengths clipped. In case
/// of standard modes like Global, Semi-Global and Local alignment, the clip operations
/// are filtered out
#[derive(Eq, PartialEq, Debug, Copy, Clone, Serialize, Deserialize)]
pub enum AlignmentOperation {
    Match,
    Subst,
    Del,
    Ins,
    Xclip(usize),
    Yclip(usize),
}

/// The modes of alignment supported by the aligner include standard modes such as
/// Global, Semi-Global and Local alignment. In addition to this, user can also invoke
/// the custom mode. In the custom mode, users can explicitly specify the clipping penalties
/// for prefix and suffix of strings 'x' and 'y' independently. Under the hood the standard
/// modes are implemented as special cases of the custom mode with the clipping penalties
/// appropriately set
#[derive(Debug, PartialEq, Eq, Copy, Clone, Serialize, Deserialize)]
pub enum AlignmentMode {
    Local,
    Semiglobal,
    Global,
    Custom,
}

/// We consider alignment between two sequences x and  y. x is the query or read sequence
/// and y is the reference or template sequence. An alignment, consisting of a score,
/// the start and end position of the alignment on sequence x and sequence y, the
/// lengths of sequences x and y, and the alignment edit operations. The start position
/// and end position of the alignment does not include the clipped regions. The length
/// of clipped regions are already encapsulated in the Alignment Operation.
#[derive(Debug, Eq, PartialEq, Clone, Serialize, Deserialize)]
pub struct Alignment {
    /// Smith-Waterman alignment score
    pub score: i32,

    /// Start position of alignment in reference
    pub ystart: usize,

    /// Start position of alignment in query
    pub xstart: usize,

    /// End position of alignment in reference
    pub yend: usize,

    /// End position of alignment in query
    pub xend: usize,

    /// Length of the reference sequence
    pub ylen: usize,

    /// Length of the query sequence
    pub xlen: usize,

    /// Vector of alignment operations
    pub operations: Vec<AlignmentOperation>,
    pub mode: AlignmentMode,
}


impl Alignment {
    /// Return the pretty formatted alignment as a String. The string
    /// contains sets of 3 lines of length 100. First line is for the
    /// sequence x, second line is for the alignment operation and the
    /// the third line is for the sequence y. A '-' in the sequence
    /// indicates a blank (insertion/deletion). The operations follow
    /// the following convention: '|' for a match, '\' for a mismatch,
    /// '+' for an insertion, 'x' for a deletion and ' ' for clipping
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alignment::pairwise::*;
    ///
    /// let x = b"CCGTCCGGCAAGGG";
    /// let y = b"AAAAACCGTTGACGGCCAA";
    /// let score = |a: u8, b: u8| if a == b {1i32} else {-2i32};
    ///
    /// let mut aligner = Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
    /// let alignment = aligner.semiglobal(x, y);
    /// println!("Semiglobal: \n{}\n", alignment.pretty(x, y));
    /// // Semiglobal:
    /// //      CCGTCCGGCAAGGG
    /// //      ||||++++\\|\||
    /// // AAAAACCGT----TGACGGCCAA

    /// let alignment = aligner.local(x, y);
    /// println!("Local: \n{}\n", alignment.pretty(x, y));
    /// // Local:
    /// //      CCGTCCGGCAAGGG
    /// //      ||||
    /// // AAAAACCGT          TGACGGCCAA

    /// let alignment = aligner.global(x, y);
    /// println!("Global: \n{}\n", alignment.pretty(x, y));
    /// // Global:
    /// // -----CCGT--CCGGCAAGGG
    /// // xxxxx||||xx\||||\|++\
    /// // AAAAACCGTTGACGGCCA--A
    /// ```
    pub fn pretty(&self, x: TextSlice, y: TextSlice) -> String {
        let mut x_pretty = String::new();
        let mut y_pretty = String::new();
        let mut inb_pretty = String::new();

        if !self.operations.is_empty() {
            let mut x_i: usize;
            let mut y_i: usize;

            // If the alignment mode is one of the standard ones, the prefix clipping is
            // implicit so we need to process it here
            match self.mode {
                AlignmentMode::Custom => {
                    x_i = 0;
                    y_i = 0;
                }
                _ => {
                    x_i = self.xstart;
                    y_i = self.ystart;
                    for k in 0..x_i {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[k]])));
                        inb_pretty.push(' ');
                        y_pretty.push(' ')
                    }
                    for k in 0..y_i {
                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[k]])));
                        inb_pretty.push(' ');
                        x_pretty.push(' ')
                    }
                }
            }

            // Process the alignment.
            for i in 0..self.operations.len() {
                match self.operations[i] {
                    AlignmentOperation::Match => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push_str("|");

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Subst => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push('\\');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Del => {
                        x_pretty.push('-');

                        inb_pretty.push('x');

                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[y_i]])));
                        y_i += 1;
                    }
                    AlignmentOperation::Ins => {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[x_i]])));
                        x_i += 1;

                        inb_pretty.push('+');

                        y_pretty.push('-');
                    }
                    AlignmentOperation::Xclip(len) => {
                        for k in 0..len {
                            x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[k]])));
                            x_i += 1;

                            inb_pretty.push(' ');

                            y_pretty.push(' ')
                        }
                    }
                    AlignmentOperation::Yclip(len) => {
                        for k in 0..len {
                            y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[k]])));
                            y_i += 1;

                            inb_pretty.push(' ');

                            x_pretty.push(' ')
                        }
                    }
                }
            }

            // If the alignment mode is one of the standard ones, the suffix clipping is
            // implicit so we need to process it here
            match self.mode {
                AlignmentMode::Custom => {}
                _ => {
                    for k in x_i..self.xlen {
                        x_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[x[k]])));
                        inb_pretty.push(' ');
                        y_pretty.push(' ')
                    }
                    for k in y_i..self.ylen {
                        y_pretty.push_str(&format!("{}", String::from_utf8_lossy(&[y[k]])));
                        inb_pretty.push(' ');
                        x_pretty.push(' ')
                    }
                }
            }
        }

        let mut s = String::new();
        let mut idx = 0;
        let step = 100; // Number of characters per line
        use std::cmp::min;

        assert_eq!(x_pretty.len(), inb_pretty.len());
        assert_eq!(y_pretty.len(), inb_pretty.len());

        let ml = x_pretty.len();

        while idx < ml {
            let rng = idx..min(idx + step, ml);
            s.push_str(&x_pretty[rng.clone()]);
            s.push_str("\n");

            s.push_str(&inb_pretty[rng.clone()]);
            s.push_str("\n");

            s.push_str(&y_pretty[rng]);
            s.push_str("\n");

            s.push_str("\n\n");
            idx += step;
        }

        s
    }

    /// Returns the optimal path in the alignment matrix
    pub fn path(&self) -> Vec<(usize, usize, AlignmentOperation)> {
        let mut path = Vec::new();

        if !self.operations.is_empty() {
            let last = match self.mode {
                AlignmentMode::Custom => (self.xlen, self.ylen),
                _ => (self.xend, self.yend),
            };
            let mut x_i = last.0;
            let mut y_i = last.1;

            let mut ops = self.operations.clone();
            ops.reverse();

            // Process the alignment.
            for i in 0..ops.len() {
                path.push((x_i, y_i, ops[i]));
                match ops[i] {
                    AlignmentOperation::Match => {
                        x_i -= 1;
                        y_i -= 1;
                    }
                    AlignmentOperation::Subst => {
                        x_i -= 1;
                        y_i -= 1;
                    }
                    AlignmentOperation::Del => {
                        y_i -= 1;
                    }
                    AlignmentOperation::Ins => {
                        x_i -= 1;
                    }
                    AlignmentOperation::Xclip(len) => {
                        x_i -= len;
                    }
                    AlignmentOperation::Yclip(len) => {
                        y_i -= len;
                    }
                }
            }
        }
        path.reverse();
        path
    }

    /// Filter out Xclip and Yclip operations from the list of operations. Useful
    /// when invoking the standard modes.
    pub fn filter_clip_operations(&mut self) {
        use self::AlignmentOperation::{Match, Subst, Ins, Del};
        self.operations
            .retain(|&ref x| (*x == Match || *x == Subst || *x == Ins || *x == Del));
    }

    /// Number of bases in reference sequence that are aligned
    pub fn y_aln_len(&self) -> usize {
        self.yend - self.ystart
    }

    /// Number of bases in query sequence that are aigned
    pub fn x_aln_len(&self) -> usize {
        self.xend - self.xstart
    }
}
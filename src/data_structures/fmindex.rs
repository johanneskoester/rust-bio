// Copyright 2014-2016 Johannes Köster, Taylor Cramer.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! The [Full-text index in Minute space index (FM-index)](https://doi.org/10.1109/SFCS.2000.892127) and
//! the FMD-Index for finding suffix array intervals matching a given pattern in linear time.
//!
//! # Examples
//!
//! ## Generate
//!
//! ```
//! use bio::alphabets::dna;
//! use bio::data_structures::bwt::{bwt, less, Occ};
//! use bio::data_structures::fmindex::{FMIndex, FMIndexable};
//! use bio::data_structures::suffix_array::suffix_array;
//!
//! let text = b"GCCTTAACATTATTACGCCTA$";
//! let alphabet = dna::n_alphabet();
//! let sa = suffix_array(text);
//! let bwt = bwt(text, &sa);
//! let less = less(&bwt, &alphabet);
//! let occ = Occ::new(&bwt, 3, &alphabet);
//! let fm = FMIndex::new(&bwt, &less, &occ);
//! ```
//!
//! ## Enclose in struct
//!
//! `FMIndex` was designed to not forcibly own the BWT and auxiliary data structures.
//! It can take a reference (`&`), owned structs or any of the more complex pointer types.
//!
//! ```
//! use bio::alphabets::dna;
//! use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
//! use bio::data_structures::fmindex::{FMIndex, FMIndexable};
//! use bio::data_structures::suffix_array::suffix_array;
//! use bio::utils::TextSlice;
//!
//! pub struct Example {
//!     fmindex: FMIndex<BWT, Less, Occ>,
//! }
//!
//! impl Example {
//!     pub fn new(text: TextSlice) -> Self {
//!         let alphabet = dna::n_alphabet();
//!         let sa = suffix_array(text);
//!         let bwt = bwt(text, &sa);
//!         let less = less(&bwt, &alphabet);
//!         let occ = Occ::new(&bwt, 3, &alphabet);
//!         let fm = FMIndex::new(bwt, less, occ);
//!         Example { fmindex: fm }
//!     }
//! }
//! ```

use std::borrow::Borrow;
use std::iter::DoubleEndedIterator;

use crate::alphabets::dna;
use crate::data_structures::bwt::{Less, Occ, BWT};
use crate::data_structures::suffix_array::SuffixArray;
use std::mem::swap;

/// A suffix array interval.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct Interval {
    pub lower: usize,
    pub upper: usize,
}

impl Interval {
    pub fn occ<SA: SuffixArray>(&self, sa: &SA) -> Vec<usize> {
        (self.lower..self.upper)
            .map(|pos| sa.get(pos).expect("Interval out of range of suffix array"))
            .collect()
    }
}

/// This enum represents the potential result states
/// from a backward_search in the fm index.  The
/// potential variants of the enum are:
/// Complete(Interval) — the query matched completely. The interval is the
/// range of suffix array indices matching the query string.
/// Partial(Intarval, usize) - some suffix of the query matched, but not the whole query.
/// The interval returned is the range of suffix array indices for the maximal
/// matching suffix, and the `usize` is the length of the maximal matching suffix.
/// Absent - None suffix of the pattern matched in the text.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub enum BackwardSearchResult {
    Complete(Interval),
    Partial(Interval, usize),
    Absent,
}

pub trait FMIndexable {
    /// Get occurrence count of symbol a in BWT[..r+1].
    fn occ(&self, r: usize, a: u8) -> usize;
    /// Also known as
    fn less(&self, a: u8) -> usize;
    fn bwt(&self) -> &BWT;

    /// Perform backward search, yielding `BackwardSearchResult` enum that
    /// contains the suffix array interval denoting exact occurrences of the given pattern
    /// of length m in the text if it exists, or the suffix array interval denoting the
    /// exact occurrences of a maximal matching suffix of the given pattern if it does
    /// not exist.  If none of the pattern can be matched, the `BackwardSearchResult` is
    /// `Absent`.
    /// Complexity: O(m).
    ///
    /// # Arguments
    ///
    /// * `pattern` - the pattern to search
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets::dna;
    /// use bio::data_structures::bwt::{bwt, less, Occ};
    /// use bio::data_structures::fmindex::{BackwardSearchResult, FMIndex, FMIndexable};
    /// use bio::data_structures::suffix_array::suffix_array;
    ///
    /// let text = b"GCCTTAACATTATTACGCCTA$";
    /// let alphabet = dna::n_alphabet();
    /// let sa = suffix_array(text);
    /// let bwt = bwt(text, &sa);
    /// let less = less(&bwt, &alphabet);
    /// let occ = Occ::new(&bwt, 3, &alphabet);
    /// let fm = FMIndex::new(&bwt, &less, &occ);
    ///
    /// let pattern = b"TTA";
    /// let bsr = fm.backward_search(pattern.iter());
    ///
    /// let positions = match bsr {
    ///     BackwardSearchResult::Complete(sai) => sai.occ(&sa),
    ///     BackwardSearchResult::Partial(sai, _l) => sai.occ(&sa),
    ///     BackwardSearchResult::Absent => Vec::<usize>::new()
    /// };
    ///
    /// assert_eq!(positions, [3, 12, 9]);
    /// ```
    fn backward_search<'b, P: Iterator<Item = &'b u8> + DoubleEndedIterator>(
        &self,
        pattern: P,
    ) -> BackwardSearchResult {
        let (mut l, mut r) = (0, self.bwt().len() - 1);
        // to keep track of the last "valid" search interval if
        // there is any valid suffix match.
        let (mut pl, mut pr) = (l, r);

        // the length of the suffix we have been able to match
        // successfully
        let mut matched_len = 0;
        // track if we exit early or not due to an empty
        // search interval.
        let mut complete_match = true;

        for &a in pattern.rev() {
            let less = self.less(a);
            pl = l;
            pr = r;
            l = less + if l > 0 { self.occ(l - 1, a) } else { 0 };
            r = less + self.occ(r, a) - 1;

            // The symbol was not found if we end up with an empty interval.
            // Terminate the LF-mapping process. In this case, also mark that
            // we do not have a complete match.
            if l > r {
                complete_match = false;
                break;
            }
            matched_len += 1;
        }

        // if we matched at least 1 character
        if matched_len > 0 {
            // if we matched the full pattern length we
            // have a complete match
            if complete_match {
                BackwardSearchResult::Complete(Interval {
                    lower: l,
                    upper: r + 1,
                })
            } else {
                // if we matched less than the full pattern length, we have
                // a partial suffix match
                BackwardSearchResult::Partial(
                    Interval {
                        lower: pl,
                        upper: pr + 1,
                    },
                    matched_len,
                )
            }
        } else {
            // if we matched nothing we have an absent result
            BackwardSearchResult::Absent
        }
    }
}

/// The Fast Index in Minute space (FM-Index, Ferragina and Manzini, 2000) for finding suffix array
/// intervals matching a given pattern.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct FMIndex<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> {
    bwt: DBWT,
    less: DLess,
    occ: DOcc,
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> FMIndexable
    for FMIndex<DBWT, DLess, DOcc>
{
    fn occ(&self, r: usize, a: u8) -> usize {
        self.occ.borrow().get(self.bwt.borrow(), r, a)
    }
    fn less(&self, a: u8) -> usize {
        self.less.borrow()[a as usize]
    }
    /// Provide a reference to the underlying BWT.
    fn bwt(&self) -> &BWT {
        self.bwt.borrow()
    }
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> FMIndex<DBWT, DLess, DOcc> {
    /// Construct a new instance of the FM index.
    ///
    /// # Arguments
    ///
    /// * `bwt` - the BWT
    /// * `less` - the less array of the BWT
    /// * `occ` - the occurence array of the BWT
    pub fn new(bwt: DBWT, less: DLess, occ: DOcc) -> Self {
        FMIndex { bwt, less, occ }
    }
}

/// A bi-interval on suffix array of the forward and reverse strand of a DNA text.
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct BiInterval {
    lower: usize,
    lower_rev: usize,
    size: usize,
    match_size: usize,
}

impl BiInterval {
    pub fn forward(&self) -> Interval {
        Interval {
            upper: self.lower + self.size,
            lower: self.lower,
        }
    }
    pub fn revcomp(&self) -> Interval {
        Interval {
            upper: self.lower_rev + self.size,
            lower: self.lower_rev,
        }
    }

    fn swapped(&self) -> BiInterval {
        BiInterval {
            lower: self.lower_rev,
            lower_rev: self.lower,
            size: self.size,
            match_size: self.match_size,
        }
    }
}

/// The FMD-Index for linear time search of supermaximal exact matches on forward and reverse
/// strand of DNA texts (Li, 2012).
#[derive(
    Default, Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize,
)]
pub struct FMDIndex<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> {
    fmindex: FMIndex<DBWT, DLess, DOcc>,
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> FMIndexable
    for FMDIndex<DBWT, DLess, DOcc>
{
    fn occ(&self, r: usize, a: u8) -> usize {
        self.fmindex.occ(r, a)
    }

    fn less(&self, a: u8) -> usize {
        self.fmindex.less(a)
    }

    /// Provide a reference to the underlying BWT.
    fn bwt(&self) -> &BWT {
        self.fmindex.bwt()
    }
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> From<FMIndex<DBWT, DLess, DOcc>>
    for FMDIndex<DBWT, DLess, DOcc>
{
    /// Construct a new instance of the FMD index (see Heng Li (2012) Bioinformatics).
    /// This expects a BWT that was created from a text over the DNA alphabet with N
    /// (`alphabets::dna::n_alphabet()`) consisting of the
    /// concatenation with its reverse complement, separated by the sentinel symbol `$`.
    /// I.e., let T be the original text and R be its reverse complement.
    /// Then, the expected text is T$R$. Further, multiple concatenated texts are allowed, e.g.
    /// T1$R1$T2$R2$T3$R3$.
    fn from(fmindex: FMIndex<DBWT, DLess, DOcc>) -> FMDIndex<DBWT, DLess, DOcc> {
        let mut alphabet = dna::n_alphabet();
        alphabet.insert(b'$');
        assert!(
            alphabet.is_word(fmindex.bwt()),
            "Expecting BWT over the DNA alphabet (including N) with the sentinel $."
        );

        FMDIndex { fmindex }
    }
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> FMDIndex<DBWT, DLess, DOcc> {
    /// Find supermaximal exact matches (of length >= l) of given pattern that overlap position i in the pattern.
    /// Complexity O(m) with pattern of length m.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets::dna;
    /// use bio::data_structures::bwt::{bwt, less, Occ};
    /// use bio::data_structures::fmindex::{FMDIndex, FMIndex};
    /// use bio::data_structures::suffix_array::suffix_array;
    ///
    /// let text = b"ATTC$GAAT$";
    /// let alphabet = dna::n_alphabet();
    /// let sa = suffix_array(text);
    /// let bwt = bwt(text, &sa);
    /// let less = less(&bwt, &alphabet);
    /// let occ = Occ::new(&bwt, 3, &alphabet);
    /// let fm = FMIndex::new(&bwt, &less, &occ);
    /// let fmdindex = FMDIndex::from(fm);
    ///
    /// let pattern = b"ATT";
    /// let intervals = fmdindex.smems(pattern, 2, 0);
    ///
    /// let forward_positions = intervals[0].0.forward().occ(&sa);
    /// let revcomp_positions = intervals[0].0.revcomp().occ(&sa);
    /// let pattern_position = intervals[0].1;
    /// let smem_len = intervals[0].2;
    ///
    /// assert_eq!(forward_positions, [0]);
    /// assert_eq!(revcomp_positions, [6]);
    /// assert_eq!(pattern_position, 0);
    /// assert_eq!(smem_len, 3);
    /// ```
    pub fn smems(&self, pattern: &[u8], i: usize, l: usize) -> Vec<(BiInterval, usize, usize)> {
        let curr = &mut Vec::new(); // pairs (biinterval, current match length)
        let prev = &mut Vec::new(); // """
        let mut matches = Vec::new(); // triples (biinterval, position on pattern, smem length)

        let mut match_len = 0;
        let mut interval = self.init_interval_with(pattern[i]);
        if interval.size != 0 {
            match_len += 1;
        }

        for &a in &pattern[i + 1..] {
            // forward extend interval
            let forward_interval = self.forward_ext(&interval, a);

            // if size changed, add last interval to list
            if interval.size != forward_interval.size {
                curr.push((interval, match_len));
            }
            // if new interval size is zero, stop, as no further forward extension is possible
            if forward_interval.size == 0 {
                break;
            }
            interval = forward_interval;
            match_len += 1;
        }
        // add the last non-zero interval
        curr.push((interval, match_len));
        // reverse intervals such that longest comes first
        curr.reverse();

        swap(curr, prev);
        let mut j = pattern.len() as isize;

        for k in (-1..i as isize).rev() {
            let a = if k == -1 { b'$' } else { pattern[k as usize] };
            curr.clear();
            // size of the last confirmed interval
            let mut last_size = -1;

            for (interval, match_len) in prev.iter() {
                // backward extend interval
                let forward_interval = self.backward_ext(interval, a);

                if (forward_interval.size == 0 || k == -1) &&
                        // interval could not be extended further
                        // if no interval has been extended this iteration,
                        // interval is maximal and can be added to the matches
                        curr.is_empty() && k < j &&
                        match_len >= &l
                {
                    j = k;
                    matches.push((*interval, (k + 1) as usize, *match_len));
                }
                // add _interval to curr (will be further extended next iteration)
                if forward_interval.size != 0 && forward_interval.size as isize != last_size {
                    last_size = forward_interval.size as isize;
                    curr.push((forward_interval, match_len + 1));
                }
            }
            if curr.is_empty() {
                break;
            }
            swap(curr, prev);
        }

        matches
    }

    /// Find all supermaximal exact matches (of length >= l) of given pattern.
    /// Complexity O(m^2) with pattern of length m.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets::dna;
    /// use bio::data_structures::bwt::{bwt, less, Occ};
    /// use bio::data_structures::fmindex::{FMDIndex, FMIndex};
    /// use bio::data_structures::suffix_array::suffix_array;
    ///
    /// let text = b"ATTCGGGG$CCCCGAAT$";
    /// let alphabet = dna::n_alphabet();
    /// let sa = suffix_array(text);
    /// let bwt = bwt(text, &sa);
    /// let less = less(&bwt, &alphabet);
    /// let occ = Occ::new(&bwt, 3, &alphabet);
    /// let fm = FMIndex::new(&bwt, &less, &occ);
    /// let fmdindex = FMDIndex::from(fm);
    ///
    /// let pattern = b"ATTGGGG";
    /// let intervals = fmdindex.all_smems(pattern, 0);
    /// assert_eq!(intervals.len(), 2);
    ///
    /// let solutions = vec![[0, 14, 0, 3], [4, 9, 3, 4]];
    /// for (i, interval) in intervals.iter().enumerate() {
    ///     let forward_positions = interval.0.forward().occ(&sa);
    ///     let revcomp_positions = interval.0.revcomp().occ(&sa);
    ///     let pattern_position = interval.1;
    ///     let smem_len = interval.2;
    ///     assert_eq!(
    ///         [
    ///             forward_positions[0],
    ///             revcomp_positions[0],
    ///             pattern_position,
    ///             smem_len
    ///         ],
    ///         solutions[i]
    ///     );
    /// }
    /// ```
    pub fn all_smems(&self, pattern: &[u8], l: usize) -> Vec<(BiInterval, usize, usize)> {
        let mut smems = Vec::new();
        let mut i0 = 0;
        while i0 < pattern.len() {
            let mut curr_smems = self.smems(pattern, i0, l);
            let mut next_i0 = i0 + 1; // this always works since:
                                      // if we have a SMEM overlapping i0, it is at least 1bp long.
                                      // If we don't have a smem, then we'll reiterate from i0+1
            for (_, p, l) in curr_smems.iter() {
                if p + l > next_i0 {
                    next_i0 = p + l;
                }
            }
            i0 = next_i0;
            smems.append(&mut curr_smems);
        }
        smems
    }

    /// Initialize interval with given start character.
    pub fn init_interval_with(&self, a: u8) -> BiInterval {
        let comp_a = dna::complement(a);
        let lower = self.fmindex.less(a);

        BiInterval {
            lower,
            lower_rev: self.fmindex.less(comp_a),
            size: self.fmindex.less(a + 1) - lower,
            match_size: 1,
        }
    }

    /// Initialize interval for empty pattern. The interval points at the whole suffix array.
    pub fn init_interval(&self) -> BiInterval {
        BiInterval {
            lower: 0,
            lower_rev: 0,
            size: self.fmindex.bwt.borrow().len(),
            match_size: 0,
        }
    }

    /// Backward extension of given interval with given character.
    pub fn backward_ext(&self, interval: &BiInterval, a: u8) -> BiInterval {
        let mut s = 0;
        let mut o = 0;
        let mut l = interval.lower_rev;
        // Interval [l(c(aP)), u(c(aP))] is a subinterval of [l(c(P)), u(c(P))] for each a,
        // starting with the lexicographically smallest ($),
        // then c(T) = A, c(G) = C, c(C) = G, N, c(A) = T, ...
        // Hence, we calculate lower revcomp bounds by iterating over
        // symbols and updating from previous one.
        for &b in b"$TGCNAtgcna" {
            l += s;
            o = if interval.lower == 0 {
                0
            } else {
                self.fmindex.occ(interval.lower - 1, b)
            };
            // calculate size
            s = self.fmindex.occ(interval.lower + interval.size - 1, b) - o;
            if b == a {
                break;
            }
        }
        // calculate lower bound
        let k = self.fmindex.less(a) + o;

        BiInterval {
            lower: k,
            lower_rev: l,
            size: s,
            match_size: interval.match_size + 1,
        }
    }

    pub fn forward_ext(&self, interval: &BiInterval, a: u8) -> BiInterval {
        let comp_a = dna::complement(a);

        self.backward_ext(&interval.swapped(), comp_a).swapped()
    }

    /// Construct a new instance of the FMD index (see Heng Li (2012) Bioinformatics)
    /// without checking whether the text is over the DNA alphabet with N.
    /// This expects a BWT that was created from a text over the DNA alphabet with N
    /// (`alphabets::dna::n_alphabet()`) consisting of the
    /// concatenation with its reverse complement, separated by the sentinel symbol `$`.
    /// I.e., let T be the original text and R be its reverse complement.
    /// Then, the expected text is T$R$. Further, multiple concatenated texts are allowed, e.g.
    /// T1$R1$T2$R2$T3$R3$.
    /// It is unsafe to construct an FMD index from an FM index that is not built on the DNA alphabet.
    pub unsafe fn from_fmindex_unchecked(
        fmindex: FMIndex<DBWT, DLess, DOcc>,
    ) -> FMDIndex<DBWT, DLess, DOcc> {
        FMDIndex { fmindex }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::alphabets::dna;
    use crate::data_structures::bwt::{bwt, less, Occ};
    use crate::data_structures::suffix_array::suffix_array;

    #[test]
    fn test_fmindex() {
        let text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);
        let fm = FMIndex::new(&bwt, &less, &occ);

        let pattern = b"TTA";
        let sai = fm.backward_search(pattern.iter());

        let positions = match sai {
            BackwardSearchResult::Complete(saint) => saint.occ(&sa),
            BackwardSearchResult::Partial(saint, _l) => saint.occ(&sa),
            BackwardSearchResult::Absent => Vec::<usize>::new(),
        };

        assert_eq!(positions, [3, 12, 9]);
    }

    #[test]
    fn test_fmindex_not_found() {
        let text = b"TCCTTAACATTATTACTCCTA$";
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);
        let fm = FMIndex::new(&bwt, &less, &occ);

        let pattern = b"TTG";
        let sai = fm.backward_search(pattern.iter());

        let positions = match sai {
            BackwardSearchResult::Complete(saint) => saint.occ(&sa),
            BackwardSearchResult::Partial(saint, _l) => saint.occ(&sa),
            BackwardSearchResult::Absent => Vec::<usize>::new(),
        };

        assert_eq!(positions, []);
    }

    #[test]
    fn test_fmindex_backward_search_optimization() {
        let text = b"GATTACA$";
        let pattern = &text[..text.len() - 1];
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);
        let fm = FMIndex::new(&bwt, &less, &occ);

        let sai = fm.backward_search(pattern.iter());

        let positions = match sai {
            BackwardSearchResult::Complete(saint) => saint.occ(&sa),
            BackwardSearchResult::Partial(saint, _l) => saint.occ(&sa),
            BackwardSearchResult::Absent => Vec::<usize>::new(),
        };

        assert_eq!(positions, [0]);
    }

    #[test]
    fn test_fmindex_backward_search_partial_match() {
        let text = b"GATTACA$";
        let pattern = b"GTACA";
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);
        let fm = FMIndex::new(&bwt, &less, &occ);

        let sai = fm.backward_search(pattern.iter());

        let mut partial_match_len = 0;
        let positions = match sai {
            BackwardSearchResult::Complete(saint) => saint.occ(&sa),
            BackwardSearchResult::Partial(saint, l) => {
                partial_match_len = l;
                saint.occ(&sa)
            }
            BackwardSearchResult::Absent => Vec::<usize>::new(),
        };

        assert_eq!(partial_match_len, 4);
        assert_eq!(positions, [3]);
    }

    #[test]
    fn test_smems() {
        let orig_text = b"GCCTTAACAT";
        let revcomp_text = dna::revcomp(orig_text);
        let text_builder: Vec<&[u8]> = vec![orig_text, b"$", &revcomp_text[..], b"$"];
        let text = text_builder.concat();

        let alphabet = dna::n_alphabet();
        let sa = suffix_array(&text);
        let bwt = bwt(&text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        let fmindex = FMIndex::new(&bwt, &less, &occ);
        let fmdindex = FMDIndex::from(fmindex);
        {
            let pattern = b"AA";
            let intervals = fmdindex.smems(pattern, 0, 0);
            let forward = intervals[0].0.forward();
            let revcomp = intervals[0].0.revcomp();
            let pattern_position = intervals[0].1;
            let smem_len = intervals[0].2;
            assert_eq!(forward.occ(&sa), [5, 16]);
            assert_eq!(revcomp.occ(&sa), [3, 14]);
            assert_eq!(pattern_position, 0);
            assert_eq!(smem_len, 2);
        }
        {
            let pattern = b"CTTAA";
            let intervals = fmdindex.smems(pattern, 1, 0);
            assert_eq!(intervals[0].0.forward().occ(&sa), [2]);
            assert_eq!(intervals[0].0.revcomp().occ(&sa), [14]);
            assert_eq!(intervals[0].1, 0);
            assert_eq!(intervals[0].2, 5);
            assert_eq!(intervals[0].0.match_size, 5);
        }
        {
            let pattern = b"CTTAA";
            let intervals = fmdindex.smems(pattern, 1, 7);
            assert!(intervals.is_empty());
        }
    }

    #[test]
    fn test_all_smems() {
        let text = b"ATTCGGGG$CCCCGAAT$";
        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);
        let fm = FMIndex::new(&bwt, &less, &occ);
        let fmdindex = FMDIndex::from(fm);

        {
            let pattern = b"ATTGGGG";
            let intervals = fmdindex.all_smems(pattern, 0);
            assert_eq!(intervals.len(), 2);
            let solutions = vec![[0, 14, 0, 3], [4, 9, 3, 4]];
            for (i, interval) in intervals.iter().enumerate() {
                let forward_positions = interval.0.forward().occ(&sa);
                let revcomp_positions = interval.0.revcomp().occ(&sa);
                let pattern_position = interval.1;
                let smem_len = interval.2;
                assert_eq!(
                    [
                        forward_positions[0],
                        revcomp_positions[0],
                        pattern_position,
                        smem_len
                    ],
                    solutions[i]
                );
            }
        }
    }

    #[test]
    fn test_init_interval() {
        let text = b"ACGT$TGCA$";

        let alphabet = dna::n_alphabet();
        let sa = suffix_array(text);
        let bwt = bwt(text, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        let fmindex = FMIndex::new(&bwt, &less, &occ);
        let fmdindex = FMDIndex::from(fmindex);
        let pattern = b"T";
        let interval = fmdindex.init_interval_with(pattern[0]);

        assert_eq!(interval.forward().occ(&sa), [3, 5]);
        assert_eq!(interval.revcomp().occ(&sa), [8, 0]);

        let empty = fmdindex.init_interval();
        let extended = fmdindex.backward_ext(&empty, pattern[0]);
        assert_eq!(extended, interval);
        let extended = fmdindex.forward_ext(&empty, pattern[0]);
        assert_eq!(extended, interval);
    }

    #[test]
    fn test_issue39() {
        let reads = b"GGCGTGGTGGCTTATGCCTGTAATCCCAGCACTTTGGGAGGTCGAAGTGGGCGG$CCGC\
                       CCACTTCGACCTCCCAAAGTGCTGGGATTACAGGCATAAGCCACCACGCC$CGAAGTGG\
                       GCGGATCACTTGAGGTCAGGAGTTGGAGACTAGCCTGGCCAACACGATGAAACCCCGTC\
                       TCTAATA$TATTAGAGACGGGGTTTCATCGTGTTGGCCAGGCTAGTCTCCAACTCCTGA\
                       CCTCAAGTGATCCGCCCACTTCG$AGCTCGAAAAATGTTTGCTTATTTTGGTAAAATTA\
                       TTCATTGACTATGCTCAGAAATCAAGCAAACTGTCCATATTTCATTTTTTG$CAAAAAA\
                       TGAAATATGGACAGTTTGCTTGATTTCTGAGCATAGTCAATGAATAATTTTACCAAAAT\
                       AAGCAAACATTTTTCGAGCT$AGCTCGAAAAATGTTTGCTTATTTTGGTAAAATTATTC\
                       ATTGACTATGCTCAGAAATCAAGCAAACTGTCCATATTTCATTTTTTGAAATTACATAT\
                       $ATATGTAATTTCAAAAAATGAAATATGGACAGTTTGCTTGATTTCTGAGCATAGTCAA\
                       TGAATAATTTTACCAAAATAAGCAAACATTTTTCGAGCT$TAAAATTTCCTCTGACAGT\
                       GTAAAAGAGATCTTCATACAAAAATCAGAATTTATATAGTCTCTTTCCAAAAGACCATA\
                       AAACCAATCAGTTAATAGTTGAT$ATCAACTATTAACTGATTGGTTTTATGGTCTTTTG\
                       GAAAGAGACTATATAAATTCTGATTTTTGTATGAAGATCTCTTTTACACTGTCAGAGGA\
                       AATTTTA$CACCTATCTACCCTGAATCTAAGTGCTAACAGGAAAGGATGCCAGATTGCA\
                       TGCCTGCTGATAAAGCCACAGTTTGGACTGTCACTCAATCACCATCGTTC$GAACGATG\
                       GTGATTGAGTGACAGTCCAAACTGTGGCTTTATCAGCAGGCATGCAATCTGGCATCCTT\
                       TCCTGTTAGCACTTAGATTCAGGGTAGATAGGTG$CATCGTTCCTCCTGTGACTCAGTA\
                       TAACAAGATTGGGAGAATACTCTACAGTTCCTGATTCCCCCACAG$CTGTGGGGGAATC\
                       AGGAACTGTAGAGTATTCTCCCAATCTTGTTATACTGAGTCACAGGAGGAACGATG$TG\
                       TAAATTCTGAGAAAAATTTGCAGGTCTTTCTTCAGGAGCATGTAATCTCTTGCTCTCTT\
                       TGTTATCTATCTATAGTACTGTAGGTTATCTGGAGTTGCT$AGCAACTCCAGATAACCT\
                       ACAGTACTATAGATAGATAACAAAGAGAGCAAGAGATTACATGCTCCTGAAGAAAGACC\
                       TGCAAATTTTTCTCAGAATTTACA$CACTTCTCCTTGTCTTTACAGACTGGTTTTGCAC\
                       TGGGAAATCCTTTCACCAGTCAGCCCAGTTAGAGATTCTG$CAGAATCTCTAACTGGGC\
                       TGACTGGTGAAAGGATTTCCCAGTGCAAAACCAGTCTGTAAAGACAAGGAGAAGTG$AA\
                       TGGAGGTATATAAATTATCTGGCAAAGTGACATATCCTGACACATTCTCCAGGATAGAT\
                       CAAATGTTAGGTCACAAAGAGAGTCTTAACAAAATT$AATTTTGTTAAGACTCTCTTTG\
                       TGACCTAACATTTGATCTATCCTGGAGAATGTGTCAGGATATGTCACTTTGCCAGATAA\
                       TTTATATACCTCCATT$TTAATTTTGTTAAGACTCTCTTTGTGACCTAACATTTGATCT\
                       ATCCTGGAGAATGTGTCAGGATATGTCACTTTGCCAGATAATTTATATACCTCCATTTT\
                       $AAAATGGAGGTATATAAATTATCTGGCAAAGTGACATATCCTGACACATTCTCCAGGA\
                       TAGATCAAATGTTAGGTCACAAAGAGAGTCTTAACAAAATTAA$TTCTTCTTTGACTCA\
                       TTGGTTGTTCAATAGTATGTTGTTTAATTTCCATATATTTGTAAATGTTTCCGTTTTCC\
                       TTCTACTATTGAATTTTTGCTTCATC$GATGAAGCAAAAATTCAATAGTAGAAGGAAAA\
                       CGGAAACATTTACAAATATATGGAAATTAAACAACATACTATTGAACAACCAATGAGTC\
                       AAAGAAGAA$AGGAAAACGGAAACATTTACAAATATATGGAAATTAAACAACATACTAT\
                       TGAACAACCAATGAGTCAAAGAAGAAATCAAAAAGAATATTAGAAAAC$GTTTTCTAAT\
                       ATTCTTTTTGATTTCTTCTTTGACTCATTGGTTGTTCAATAGTATGTTGTTTAATTTCC\
                       ATATATTTGTAAATGTTTCCGTTTTCCT$TTAGAAAACAAGCTGACAAAAAAATAAAAA\
                       AACACAACATAGCAAAACTTAGAAATGCAGCAAAGGCAGTACTAAAGAGGGAAATTTAT\
                       AGCAATAAATGC$GCATTTATTGCTATAAATTTCCCTCTTTAGTACTGCCTTTGCTGCA\
                       TTTCTAAGTTTTGCTATGTTGTGTTTTTTTATTTTTTTGTCAGCTTGTTTTCTAA$TTT\
                       ATTGCTATAAATTTCCCTCTTTAGTACTGCCTTTGCTGCATTTCTAAGTTTTGCTATGT\
                       TGTGTTTTTTTATTTTTTTGTCAGCTTGTTTTCTA$TAGAAAACAAGCTGACAAAAAAA\
                       TAAAAAAACACAACATAGCAAAACTTAGAAATGCAGCAAAGGCAGTACTAAAGAGGGAA\
                       ATTTATAGCAATAAA$TCTTTCTTCTTTTTTAAGGTAGGCATTTATTGCTATAAATTTC\
                       CCTCTTTAGTACTGCCTTTG$CAAAGGCAGTACTAAAGAGGGAAATTTATAGCAATAAA\
                       TGCCTACCTTAAAAAAGAAGAAAGA$";

        let alphabet = dna::n_alphabet();
        let sa = suffix_array(reads);
        let bwt = bwt(reads, &sa);
        let less = less(&bwt, &alphabet);
        let occ = Occ::new(&bwt, 3, &alphabet);

        let fmindex = FMIndex::new(&bwt, &less, &occ);
        let fmdindex = FMDIndex::from(fmindex);

        let read = b"GGCGTGGTGGCTTATGCCTGTAATCCCAGCACTTTGGGAGGTCGAAGTGGGCGG";
        let read_pos = 0;

        for i in 0..read.len() {
            println!("i {}", i);
            let intervals = fmdindex.smems(read, i, 0);
            println!("{:?}", intervals);
            let matches = intervals
                .iter()
                .flat_map(|interval| interval.0.forward().occ(&sa))
                .collect::<Vec<usize>>();
            assert_eq!(matches, vec![read_pos]);
        }
    }
}

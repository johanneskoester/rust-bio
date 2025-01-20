// Copyright 2014-2016 Johannes Köster, Taylor Cramer.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Suffix arrays and related algorithms.
//! The implementation is based on the lecture notes
//! "Algorithmen auf Sequenzen", Kopczynski, Marschall, Martin and Rahmann, 2008 - 2015.
//! The original algorithm desciption can be found in:
//! [Ge Nong, Sen Zhang, Wai Hong Chan: Two Efficient Algorithms for Linear Time Suffix Array Construction. IEEE Trans. Computers 60(10): 1471–1484 (2011)](https://doi.org/10.1109/TC.2010.188)
//!
//! # Examples
//!
//! ```
//! use bio::data_structures::suffix_array::suffix_array;
//! let text = b"GCCTTAACATTATTACGCCTA$";
//! let pos = suffix_array(text);
//! assert_eq!(
//!     pos,
//!     vec![21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18, 2, 16, 0, 19, 4, 13, 10, 3, 12, 9]
//! );
//! ```

use std::borrow::Borrow;
use std::cmp;
use std::collections::HashMap;
use std::fmt::Debug;
use std::hash::BuildHasherDefault;
use std::iter;
use std::ops::Deref;

use num_integer::Integer;
use num_traits::{cast, NumCast, Unsigned};

use bv::{BitVec, Bits, BitsMut};
use vec_map::VecMap;

use fxhash::FxHasher;

use crate::alphabets::{Alphabet, RankTransform};
use crate::data_structures::bwt::{Less, Occ, BWT};
use crate::data_structures::smallints::SmallInts;

pub type LCPArray = SmallInts<i8, isize>;
pub type RawSuffixArray = Vec<usize>;
pub type RawSuffixArraySlice<'a> = &'a [usize];

type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// A trait exposing general functionality of suffix arrays.
pub trait SuffixArray {
    fn get(&self, index: usize) -> Option<usize>;
    fn len(&self) -> usize;
    fn is_empty(&self) -> bool;

    /// Sample the suffix array with the given sample rate.
    ///
    /// # Arguments
    ///
    /// * `text` - text that the suffix array is built on
    /// * `bwt` - the corresponding BWT
    /// * `less` - the corresponding less array
    /// * `occ` - the corresponding occ table
    /// * `sampling_rate` - if sampling rate is k, every k-th entry will be kept
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets::dna;
    /// use bio::data_structures::bwt::{bwt, less, Occ};
    /// use bio::data_structures::suffix_array::{suffix_array, SuffixArray};
    ///
    /// let text = b"ACGCGAT$";
    /// let alphabet = dna::n_alphabet();
    /// let sa = suffix_array(text);
    /// let bwt = bwt(text, &sa);
    /// let less = less(&bwt, &alphabet);
    /// let occ = Occ::new(&bwt, 3, &alphabet);
    /// let sampled = sa.sample(text, &bwt, &less, &occ, 2);
    ///
    /// for i in 0..sa.len() {
    ///     assert_eq!(sa.get(i), sampled.get(i));
    /// }
    /// ```
    fn sample<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>>(
        &self,
        text: &[u8],
        bwt: DBWT,
        less: DLess,
        occ: DOcc,
        sampling_rate: usize,
    ) -> SampledSuffixArray<DBWT, DLess, DOcc> {
        let mut sample =
            Vec::with_capacity((self.len() as f32 / sampling_rate as f32).ceil() as usize);
        let mut extra_rows = HashMapFx::default();
        let sentinel = sentinel(text);

        for i in 0..self.len() {
            let idx = self.get(i).unwrap();
            if (i % sampling_rate) == 0 {
                sample.push(idx);
            } else if bwt.borrow()[i] == sentinel {
                // If bwt lookup will return a sentinel
                // Text suffixes that begin right after a sentinel are always saved as extra rows
                // to help deal with FM index last to front inaccuracy when there are many sentinels
                extra_rows.insert(i, idx);
            }
        }

        SampledSuffixArray {
            bwt,
            less,
            occ,
            sample,
            s: sampling_rate,
            extra_rows,
            sentinel,
        }
    }
}

/// A sampled suffix array.
#[derive(Default, Clone, Eq, PartialEq, Debug, Serialize, Deserialize)]
pub struct SampledSuffixArray<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> {
    bwt: DBWT,
    less: DLess,
    occ: DOcc,
    sample: Vec<usize>,
    s: usize, // Rate of sampling
    extra_rows: HashMapFx<usize, usize>,
    sentinel: u8,
}

impl SuffixArray for RawSuffixArray {
    fn get(&self, index: usize) -> Option<usize> {
        // Explicitly written out because Vec::get(index) generates a recursion warning
        if index < self.len() {
            Some(self[index])
        } else {
            None
        }
    }

    fn len(&self) -> usize {
        Vec::len(self)
    }

    fn is_empty(&self) -> bool {
        Vec::is_empty(self)
    }
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>> SuffixArray
    for SampledSuffixArray<DBWT, DLess, DOcc>
{
    fn get(&self, index: usize) -> Option<usize> {
        if index < self.len() {
            let mut pos = index;
            let mut offset = 0;
            loop {
                if pos % self.s == 0 {
                    return Some(self.sample[pos / self.s] + offset);
                }

                let c = self.bwt.borrow()[pos];

                if c == self.sentinel {
                    // Check if next character in the bwt is the sentinel
                    // If so, there must be a cached result to workaround FM index last to front
                    // mapping inaccuracy when there are multiple sentinels
                    // This branch should rarely be triggered so the performance impact
                    // of hashmap lookups would be low
                    return Some(self.extra_rows[&pos] + offset);
                }

                pos = self.less.borrow()[c as usize]
                    + self.occ.borrow().get(self.bwt.borrow(), pos - 1, c);
                offset += 1;
            }
        } else {
            None
        }
    }

    fn len(&self) -> usize {
        self.bwt.borrow().len()
    }

    fn is_empty(&self) -> bool {
        self.bwt.borrow().is_empty()
    }
}

impl<DBWT: Borrow<BWT>, DLess: Borrow<Less>, DOcc: Borrow<Occ>>
    SampledSuffixArray<DBWT, DLess, DOcc>
{
    /// Get the sampling rate of the suffix array.
    pub fn sampling_rate(&self) -> usize {
        self.s
    }

    /// Get a reference to the internal BWT.
    pub fn bwt(&self) -> &BWT {
        self.bwt.borrow()
    }

    /// Get a reference to the internal Less.
    pub fn less(&self) -> &Less {
        self.less.borrow()
    }

    /// Get a reference to the internal Occ.
    pub fn occ(&self) -> &Occ {
        self.occ.borrow()
    }
}

/// Construct suffix array for given text of length n.
/// Complexity: O(n).
/// This is an implementation of the induced sorting as presented by
/// Ge Nong, Sen Zhang und Wai Hong Chan (2009), also known as SAIS.
/// The implementation is based on the following lecture notes:
/// http://ls11-www.cs.tu-dortmund.de/people/rahmann/algoseq.pdf
///
/// The idea is to first mark positions as L or S, with L being a position
/// the suffix of which is lexicographically larger than that of the next position.
/// Then, LMS-positions (leftmost S) are S-positions right to an L-position.
/// An LMS substring is the substring from one LMS position to the next (inclusive).
/// The algorithm works as follows:
///
/// 1. Sort LMS positions: first step 2 is applied to the unsorted sequence
///    of positions. Surprisingly, this sorts the LMS substrings. If all substrings
///    are different, LMS positions (and their suffixes) are sorted. Else, a reduced
///    text is build (at most half the size of the original text) and we recurse into
///    suffix array construction on the reduced text, yielding the sorted LMS positions.
/// 2. Derive the order of the other positions/suffixes from the (sorted) LMS positions.
///    For this, the (still empty) suffix array is partitioned into buckets.
///    Each bucket denotes an interval of suffixes with the same first symbol.
///    We know that the L-suffixes have to occur first in the buckets, because they
///    have to be lexicographically smaller than the S-suffixes with the same first letter.
///    The LMS-positions can now be used to insert the L-positions in the correct order
///    into the buckets.
///    Then, the S-positions can be inserted, again using the already existing entries
///    in the array.
///
/// # Arguments
///
/// * `text` - the text, ended by sentinel symbol (being lexicographically smallest). The text may
///   also contain multiple sentinel symbols, used to concatenate multiple sequences without mixing
///   their suffixes together.
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::suffix_array;
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
/// assert_eq!(
///     pos,
///     vec![21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18, 2, 16, 0, 19, 4, 13, 10, 3, 12, 9]
/// );
/// ```
pub fn suffix_array(text: &[u8]) -> RawSuffixArray {
    let n = text.len();
    let alphabet = Alphabet::new(text);
    let sentinel_count = sentinel_count(text);
    let mut sais = Sais::new(n);

    match alphabet.len() + sentinel_count {
        a if a <= std::u8::MAX as usize => {
            sais.construct(&transform_text::<u8>(text, &alphabet, sentinel_count))
        }
        a if a <= std::u16::MAX as usize => {
            sais.construct(&transform_text::<u16>(text, &alphabet, sentinel_count))
        }
        a if a <= std::u32::MAX as usize => {
            sais.construct(&transform_text::<u32>(text, &alphabet, sentinel_count))
        }
        _ => sais.construct(&transform_text::<u64>(text, &alphabet, sentinel_count)),
    }

    sais.pos
}

/// Construct suffix array for given text from integer alphabet.
/// Complexity: O(n).
/// # Arguments
///
/// * `text` - the text, ended by sentinel symbol (being lexicographically smallest).
/// All symbols, from lexicographically smallest to largest, need to be present in the text,
/// otherwise SAIS algorithm panics on 'index out of bounds' error.
/// The text may also contain multiple sentinel symbols, used to concatenate
/// multiple sequences without mixing their suffixes together.
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::suffix_array_int;
/// let text: Vec<usize> = vec![3, 2, 2, 4, 4, 1, 2, 1, 0];
/// let sa = suffix_array_int(&text);
/// assert_eq!(sa, vec![8, 7, 5, 6, 1, 2, 0, 4, 3]);
/// ```
pub fn suffix_array_int<T>(text: &[T]) -> RawSuffixArray
where
    T: Integer + Unsigned + NumCast + Copy + Debug,
{
    let mut sais = Sais::new(text.len());
    sais.construct(&text);
    sais.pos
}

/// Construct lcp array for given text and suffix array of length n.
/// Complexity: O(n).
///
/// # Arguments
///
/// * `text` - the text ended by sentinel symbol (being lexicographically smallest)
/// * `pos` - the suffix array for the text
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::{lcp, suffix_array};
/// let text = b"GCCTTAACATTATTACGCCTA$";
/// let pos = suffix_array(text);
///
/// // obtain compressed LCP array
/// let lcp = lcp(text, &pos);
///
/// // get most values in O(1).
/// assert_eq!(lcp.get(6).unwrap(), 4);
///
/// // obtain uncompressed LCP array.
/// let uncompressed = lcp.decompress();
/// assert_eq!(
///     uncompressed,
///     [-1, 0, 1, 1, 2, 1, 4, 0, 1, 3, 1, 1, 2, 0, 4, 0, 2, 2, 2, 1, 3, 3, -1]
/// )
/// ```
pub fn lcp<SA: Deref<Target = RawSuffixArray>>(text: &[u8], pos: SA) -> LCPArray {
    assert_eq!(text.len(), pos.len());
    let n = text.len();

    // provide the lexicographical rank for each suffix
    let mut rank: Vec<usize> = iter::repeat(0).take(n).collect();
    for (r, p) in pos.iter().enumerate() {
        rank[*p] = r;
    }

    let mut lcp = SmallInts::from_elem(-1, n + 1);
    let mut l = 0usize;
    for (p, &r) in rank.iter().enumerate().take(n - 1) {
        // since the sentinel has rank 0 and is excluded above,
        // we will never have a negative index below
        let pred = pos[r - 1];
        while pred + l < n && p + l < n && text[p + l] == text[pred + l] {
            l += 1;
        }
        lcp.set(r, l as isize);
        l = if l > 0 { l - 1 } else { 0 };
    }

    lcp
}

/// Calculate all locally shortest unique substrings from a given suffix and lcp array
/// (Ohlebusch (2013). "Bioinformatics Algorithms". ISBN 978-3-00-041316-2).
/// Complexity: O(n)
///
/// # Arguments
///
/// * `pos` - the suffix array
/// * `lcp` - the lcp array
///
/// # Returns
///
/// An vector of the length of the shortest unique substring for each position of the text.
/// Suffixes are excluded. If no unique substring starts at a given position, the entry is `None`.
///
/// # Example
///
/// ```
/// use bio::data_structures::suffix_array::{lcp, shortest_unique_substrings, suffix_array};
/// let text = b"GCTGCTA$";
/// let pos = suffix_array(text);
///
/// // obtain compressed LCP array
/// let lcp = lcp(text, &pos);
///
/// // calculate shortest unique substrings
/// let sus = shortest_unique_substrings(&pos, &lcp);
/// assert_eq!(
///     sus,
///     [
///         Some(4),
///         Some(3),
///         Some(2),
///         Some(4),
///         Some(3),
///         Some(2),
///         Some(1),
///         Some(1)
///     ]
/// );
/// ```
pub fn shortest_unique_substrings<SA: SuffixArray>(pos: &SA, lcp: &LCPArray) -> Vec<Option<usize>> {
    let n = pos.len();
    // Initialize array representing the length of the shortest unique substring starting at position i
    let mut sus = vec![None; n];
    for i in 0..n {
        // The longest common prefixes (LCP) of suffix pos[i] with its predecessor and successor are not unique.
        // In turn the their maximum + 1 is the length of the shortest unique substring starting at pos[i].
        let len = 1 + cmp::max(lcp.get(i).unwrap(), lcp.get(i + 1).unwrap_or(0)) as usize;
        let p = pos.get(i).unwrap();
        // Check if the suffix pos[i] is a prefix of pos[i+1]. In that case, there is no unique substring
        // at this position.
        if n - p >= len {
            sus[p] = Some(len);
        }
    }
    sus
}

/// Return last character of the text (expected to be the sentinel).
fn sentinel(text: &[u8]) -> u8 {
    text[text.len() - 1]
}

/// Count the sentinels occurring in the text given that the last character is the sentinel.
fn sentinel_count(text: &[u8]) -> usize {
    let sentinel = sentinel(text);
    assert!(
        text.iter().all(|&a| a >= sentinel),
        "Expecting extra sentinel symbol being lexicographically smallest at the end of the \
         text."
    );

    text.iter()
        .fold(0, |count, &a| count + (a == sentinel) as usize)
}

/// Transform the given text into integers for usage in `SAIS`.
fn transform_text<T: Integer + Unsigned + NumCast + Copy + Debug>(
    text: &[u8],
    alphabet: &Alphabet,
    sentinel_count: usize,
) -> Vec<T> {
    let sentinel = sentinel(text);
    let transform = RankTransform::new(alphabet);
    let offset = sentinel_count - 1;

    let mut transformed: Vec<T> = Vec::with_capacity(text.len());
    let mut s = sentinel_count;
    for &a in text.iter() {
        if a == sentinel {
            s -= 1;
            transformed.push(cast(s).unwrap());
        } else {
            transformed
                .push(cast(*(transform.ranks.get(a as usize)).unwrap() as usize + offset).unwrap());
        }
    }

    transformed
}

/// SAIS implementation (see function `suffix_array` for description).
struct Sais {
    pos: Vec<usize>,
    lms_pos: Vec<usize>,
    reduced_text_pos: Vec<usize>,
    bucket_sizes: VecMap<usize>,
    bucket_start: Vec<usize>,
    bucket_end: Vec<usize>,
}

impl Sais {
    /// Create a new instance.
    fn new(n: usize) -> Self {
        Sais {
            pos: Vec::with_capacity(n),
            lms_pos: Vec::with_capacity(n),
            reduced_text_pos: vec![0; n],
            bucket_sizes: VecMap::new(),
            bucket_start: Vec::with_capacity(n),
            bucket_end: Vec::with_capacity(n),
        }
    }

    /// Init buckets.
    fn init_bucket_start<T: Integer + Unsigned + NumCast + Copy>(&mut self, text: &[T]) {
        self.bucket_sizes.clear();
        self.bucket_start.clear();

        for &c in text.iter() {
            if !self.bucket_sizes.contains_key(cast(c).unwrap()) {
                self.bucket_sizes.insert(cast(c).unwrap(), 0);
            }
            *(self.bucket_sizes.get_mut(cast(c).unwrap()).unwrap()) += 1;
        }

        let mut sum = 0;
        for &size in self.bucket_sizes.values() {
            self.bucket_start.push(sum);
            sum += size;
        }
    }

    /// Initialize pointers to the last element of the buckets.
    fn init_bucket_end<T: Integer + Unsigned + NumCast + Copy>(&mut self, text: &[T]) {
        self.bucket_end.clear();
        for &r in self.bucket_start[1..].iter() {
            self.bucket_end.push(r - 1);
        }
        self.bucket_end.push(text.len() - 1);
    }

    /// Check if two LMS substrings are equal.
    fn lms_substring_eq<T: Integer + Unsigned + NumCast + Copy>(
        &self,
        text: &[T],
        pos_types: &PosTypes,
        i: usize,
        j: usize,
    ) -> bool {
        for k in 0.. {
            let lmsi = pos_types.is_lms_pos(i + k);
            let lmsj = pos_types.is_lms_pos(j + k);
            if text[i + k] != text[j + k] {
                // different symbols
                return false;
            }
            if lmsi != lmsj {
                // different length
                return false;
            }
            if k > 0 && lmsi && lmsj {
                // same symbols and same length
                return true;
            }
        }
        false
    }

    /// Sort LMS suffixes.
    fn sort_lms_suffixes<
        T: Integer + Unsigned + NumCast + Copy + Debug,
        S: Integer + Unsigned + NumCast + Copy + Debug,
    >(
        &mut self,
        text: &[T],
        pos_types: &PosTypes,
        lms_substring_count: usize,
    ) {
        // if less than 2 LMS substrings are present, no further sorting is needed
        if lms_substring_count > 1 {
            // sort LMS suffixes by recursively building SA on reduced text
            let mut reduced_text: Vec<S> = vec![cast(0).unwrap(); lms_substring_count];
            let mut label = 0;
            reduced_text[self.reduced_text_pos[self.pos[0]]] = cast(label).unwrap();
            let mut prev = None;
            for &p in &self.pos {
                if pos_types.is_lms_pos(p) {
                    // choose same label if substrings are equal
                    if prev.is_some() && !self.lms_substring_eq(text, pos_types, prev.unwrap(), p) {
                        label += 1;
                    }
                    reduced_text[self.reduced_text_pos[p]] = cast(label).unwrap();
                    prev = Some(p);
                }
            }

            // if we have less labels than substrings, we have to sort by recursion
            // because two or more substrings are equal
            if label + 1 < lms_substring_count {
                // backup lms_pos
                let lms_pos = self.lms_pos.clone();
                // recurse SA construction for reduced text
                self.construct(&reduced_text);
                // obtain sorted lms suffixes
                self.lms_pos.clear();
                for &p in &self.pos {
                    self.lms_pos.push(lms_pos[p]);
                }
            } else {
                // otherwise, lms_pos is updated with the sorted suffixes from pos
                // obtain sorted lms suffixes
                self.lms_pos.clear();
                for &p in &self.pos {
                    if pos_types.is_lms_pos(p) {
                        self.lms_pos.push(p);
                    }
                }
            }
        }
    }

    /// Construct the suffix array.
    fn construct<T: Integer + Unsigned + NumCast + Copy + Debug>(&mut self, text: &[T]) {
        let pos_types = PosTypes::new(text);
        self.calc_lms_pos(text, &pos_types);
        self.calc_pos(text, &pos_types);
    }

    /// Step 1 of the SAIS algorithm.
    fn calc_lms_pos<T: Integer + Unsigned + NumCast + Copy + Debug>(
        &mut self,
        text: &[T],
        pos_types: &PosTypes,
    ) {
        let n = text.len();

        // collect LMS positions
        self.lms_pos.clear();
        let mut i = 0;
        for r in 0..n {
            if pos_types.is_lms_pos(r) {
                self.lms_pos.push(r);
                self.reduced_text_pos[r] = i;
                i += 1;
            }
        }

        // sort LMS substrings by applying step 2 with unsorted LMS positions
        self.calc_pos(text, pos_types);

        let lms_substring_count = self.lms_pos.len();

        if lms_substring_count <= std::u8::MAX as usize {
            self.sort_lms_suffixes::<T, u8>(text, pos_types, lms_substring_count);
        } else if lms_substring_count <= std::u16::MAX as usize {
            self.sort_lms_suffixes::<T, u16>(text, pos_types, lms_substring_count);
        } else if lms_substring_count <= std::u32::MAX as usize {
            self.sort_lms_suffixes::<T, u32>(text, pos_types, lms_substring_count);
        } else {
            self.sort_lms_suffixes::<T, u64>(text, pos_types, lms_substring_count);
        }
    }

    /// Step 2 of the SAIS algorithm.
    fn calc_pos<T: Integer + Unsigned + NumCast + Copy>(
        &mut self,
        text: &[T],
        pos_types: &PosTypes,
    ) {
        let n = text.len();
        self.pos.clear();

        self.init_bucket_start(text);
        self.init_bucket_end(text);

        // init all positions as unknown (n-1 is max position)
        self.pos.resize(n, n);

        // insert LMS positions to the end of their buckets
        for &p in self.lms_pos.iter().rev() {
            let c: usize = cast(text[p]).unwrap();
            self.pos[self.bucket_end[c]] = p;
            // subtract without overflow: last -1 will cause overflow, but it does not matter
            self.bucket_end[c] = self.bucket_end[c].wrapping_sub(1);
        }

        // reset bucket ends
        self.init_bucket_end(text);

        // insert L-positions into buckets
        for r in 0..n {
            let p = self.pos[r];
            // ignore undefined positions and the zero since it has no predecessor
            if p == n || p == 0 {
                continue;
            }
            let pred = p - 1;
            if pos_types.is_l_pos(pred) {
                let c: usize = cast(text[pred]).unwrap();
                self.pos[self.bucket_start[c]] = pred;
                self.bucket_start[c] += 1;
            }
        }

        // insert S-positions into buckets
        for r in (0..n).rev() {
            let p = self.pos[r];
            if p == 0 {
                continue;
            }
            let pred = p - 1;
            if pos_types.is_s_pos(pred) {
                let c: usize = cast(text[pred]).unwrap();
                self.pos[self.bucket_end[c]] = pred;
                // subtract without overflow: last -1 will cause overflow, but it won't be used
                self.bucket_end[c] = self.bucket_end[c].wrapping_sub(1);
            }
        }
    }
}

/// Position types (L or S).
#[derive(Debug)]
struct PosTypes {
    pos_types: BitVec,
}

impl PosTypes {
    /// Calculate the text position type.
    /// L-type marks suffixes being lexicographically larger than their successor,
    /// S-type marks the others.
    /// This function fills a BitVec, with 1-bits denoting S-type
    /// and 0-bits denoting L-type.
    ///
    /// # Arguments
    ///
    /// * `text` - the text, ending with a sentinel.
    fn new<T: Integer + Unsigned + NumCast + Copy>(text: &[T]) -> Self {
        let n = text.len();
        let mut pos_types = BitVec::new_fill(false, n as u64);
        pos_types.set_bit(n as u64 - 1, true);

        for p in (0..n - 1).rev() {
            if text[p] == text[p + 1] {
                // if the characters are equal, the next position determines
                // the lexicographical order
                let v = pos_types.get_bit(p as u64 + 1);
                pos_types.set_bit(p as u64, v);
            } else {
                pos_types.set_bit(p as u64, text[p] < text[p + 1]);
            }
        }

        PosTypes { pos_types }
    }

    /// Check if p is S-position.
    fn is_s_pos(&self, p: usize) -> bool {
        self.pos_types.get_bit(p as u64)
    }

    /// Check if p is L-position.
    fn is_l_pos(&self, p: usize) -> bool {
        !self.pos_types.get_bit(p as u64)
    }

    /// Check if p is LMS-position.
    fn is_lms_pos(&self, p: usize) -> bool {
        p != 0 && self.is_s_pos(p) && self.is_l_pos(p - 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::{transform_text, PosTypes, Sais};
    use crate::alphabets::{dna, Alphabet};
    use crate::data_structures::bwt::{bwt, less};
    use bv::{BitVec, BitsPush};
    use rand;
    use rand::prelude::*;
    use std::str;

    #[test]
    fn test_pos_types() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let pos_types = PosTypes::new(&text);
        //let mut test = BitSlice::from_slice(&[0b01100110, 0b10010011, 0b01100100]).to_owned();
        let mut test = BitVec::new();
        test.push_block(0b001001101100100101100110);
        test.truncate(n as u64);
        assert_eq!(pos_types.pos_types, test);
        let lms_pos: Vec<usize> = (0..n).filter(|&p| pos_types.is_lms_pos(p)).collect();
        assert_eq!(lms_pos, vec![1, 5, 8, 11, 14, 17, 21]);
    }

    #[test]
    fn test_buckets() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let mut sais = Sais::new(n);
        sais.init_bucket_start(&text);
        assert_eq!(sais.bucket_start, vec![0, 1, 7, 13, 15]);
        sais.init_bucket_end(&text);
        assert_eq!(sais.bucket_end, vec![0, 6, 12, 14, 21]);
    }

    #[test]
    fn test_pos() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let mut sais = Sais::new(n);
        let pos_types = PosTypes::new(&text);
        sais.lms_pos = vec![21, 5, 14, 8, 11, 17, 1];
        sais.calc_pos(&text, &pos_types);
        assert_eq!(
            sais.pos,
            vec![21, 20, 5, 6, 14, 11, 8, 7, 17, 1, 15, 18, 2, 16, 0, 19, 4, 13, 10, 3, 12, 9,]
        );
    }

    #[test]
    fn test_lms_pos() {
        let orig_text = b"GCCTTAACATTATTACGCCTA$";
        let alphabet = Alphabet::new(orig_text);
        let text: Vec<u8> = transform_text(orig_text, &alphabet, 1);
        let n = text.len();

        let mut sais = Sais::new(n);
        let pos_types = PosTypes::new(&text);
        sais.calc_lms_pos(&text, &pos_types);
    }

    #[test]
    fn test_issue10_1() {
        let text = b"TGTGTGTGTG$";
        let pos = suffix_array(text);
        assert_eq!(pos, [10, 9, 7, 5, 3, 1, 8, 6, 4, 2, 0]);
    }

    #[test]
    fn test_issue10_2() {
        let text = b"TGTGTGTG$";
        let pos = suffix_array(text);
        assert_eq!(pos, [8, 7, 5, 3, 1, 6, 4, 2, 0]);
    }

    #[test]
    fn test_handles_sentinels_properly() {
        let reads = b"TACTCCGCTAGGGACACCTAAATAGATACTCGCAAAGGCGACTGATATATCCTTAGGTCGAAGAGATACCAGAGAAATAGTAGGTCTTAGGCTAGTCCTT$AAGGACTAGCCTAAGACCTACTATTTCTCTGGTATCTCTTCGACCTAAGGATATATCAGTCGCCTTTGCGAGTATCTATTTAGGTGTCCCTAGCGGAGTA$TAGGGACACCTAAATAGATACTCGCAAAGGCGACTGATATATCCTTAGGTCGAAGAGATACCAGAGAAATAGTAGGTCTTAGGCTAGTCCTTGTCCAGTA$TACTGGACAAGGACTAGCCTAAGACCTACTATTTCTCTGGTATCTCTTCGACCTAAGGATATATCAGTCGCCTTTGCGAGTATCTATTTAGGTGTCCCTA$ACGCACCCCGGCATTCGTCGACTCTACACTTAGTGGAACATACAAATTCGCTCGCAGGAGCGCCTCATACATTCTAACGCAGTGATCTTCGGCTGAGACT$AGTCTCAGCCGAAGATCACTGCGTTAGAATGTATGAGGCGCTCCTGCGAGCGAATTTGTATGTTCCACTAAGTGTAGAGTCGACGAATGCCGGGGTGCGT$";
        suffix_array(reads);
    }

    fn str_from_pos(sa: &[usize], text: &[u8], index: usize) -> String {
        String::from(
            str::from_utf8(&text[sa[index]..])
                .unwrap()
                .split('$')
                .next()
                .unwrap_or(""),
        ) + "$"
    }

    fn rand_seqs(num_seqs: usize, seq_len: usize) -> Vec<u8> {
        let mut rng = rand::thread_rng();
        let alpha = [b'A', b'T', b'C', b'G', b'N'];
        let seqs = (0..num_seqs)
            .map(|_| {
                let len = rng.gen_range((seq_len / 2)..=seq_len);
                (0..len)
                    .map(|_| *alpha.choose(&mut rng).unwrap())
                    .collect::<Vec<u8>>()
            })
            .collect::<Vec<_>>();
        let mut res = seqs.join(&b'$');
        res.push(b'$');
        res
    }

    #[test]
    fn test_sorts_lexically() {
        let mut test_cases = vec![(&b"A$C$G$T$"[..], "simple"),
             (&b"A$A$T$T$"[..], "duplicates"),
             (&b"AA$GA$CA$TA$TC$TG$GT$GC$"[..], "two letter"),
             (&b"AGCCAT$CAGCC$"[..], "substring"),
             (&b"GTAGGCCTAATTATAATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAA$\
                AATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAATGGCTATTCCAATA$"[..],
                "complex"),
             (&b"GTAGGCCTAATTATAATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAA$\
                TTCGACGCTGACCTCTTGAGGTTCCATTACCCGGCTACTGATGCTAAAATCCTGGCAGCCCGAGCAATACGAAATGTCCGCTGATTATAATTAGGCCTAC$\
                AATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAATGGCTATTCCAATA$\
                TATTGGAATAGCCATTCGACGCTGACCTCTTGAGGTTCCATTACCCGGCTACTGATGCTAAAATCCTGGCAGCCCGAGCAATACGAAATGTCCGCTGATT$"[..],
                "complex with revcomps"),
             ];
        let num_rand = 100;
        let rand_cases: Vec<_> = (0..num_rand).map(|i| rand_seqs(10, i * 10)).collect();
        test_cases.extend(
            rand_cases
                .iter()
                .map(|case| (case.as_ref(), "rand test case")),
        );

        for &(text, test_name) in test_cases.iter() {
            let pos = suffix_array(text);
            for i in 0..(pos.len() - 2) {
                // Check that every element in the suffix array is lexically <= the next elem
                let cur = str_from_pos(&pos, text, i);
                let next = str_from_pos(&pos, text, i + 1);

                assert!(
                    cur <= next,
                    "Failed:\n{}\n{}\nat positions {} and {} are out of order in \
                         test: {}",
                    cur,
                    next,
                    pos[i],
                    pos[i + 1],
                    test_name
                );
            }
        }
    }

    #[test]
    fn test_sampled_matches() {
        let mut test_cases = vec![(&b"A$C$G$T$"[..], "simple"),
             (&b"A$A$T$T$"[..], "duplicates"),
             (&b"AA$GA$CA$TA$TC$TG$GT$GC$"[..], "two letter"),
             (&b"AGCCAT$\
                CAGCC$"[..],
                "substring"),
             (&b"GTAGGCCTAATTATAATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAA$\
                AATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAATGGCTATTCCAATA$"[..],
                "complex"),
             (&b"GTAGGCCTAATTATAATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAA$\
                TTCGACGCTGACCTCTTGAGGTTCCATTACCCGGCTACTGATGCTAAAATCCTGGCAGCCCGAGCAATACGAAATGTCCGCTGATTATAATTAGGCCTAC$\
                AATCAGCGGACATTTCGTATTGCTCGGGCTGCCAGGATTTTAGCATCAGTAGCCGGGTAATGGAACCTCAAGAGGTCAGCGTCGAATGGCTATTCCAATA$\
                TATTGGAATAGCCATTCGACGCTGACCTCTTGAGGTTCCATTACCCGGCTACTGATGCTAAAATCCTGGCAGCCCGAGCAATACGAAATGTCCGCTGATT$"[..],
                "complex with revcomps"),
             (&b"GTAG$GCCTAAT$TATAATCAG$"[..], "issue70"),
             (&b"TGTGTGTGTG$"[..], "repeating"),
             (&b"TACTCCGCTAGGGACACCTAAATAGATACTCGCAAAGGCGACTGATATATCCTTAGGTCGAAGAGATACCAGAGAAATAGTAGGTCTTAGGCTAGTCCTT$AAGGACTAGCCTAAGACCTACTATTTCTCTGGTATCTCTTCGACCTAAGGATATATCAGTCGCCTTTGCGAGTATCTATTTAGGTGTCCCTAGCGGAGTA$TAGGGACACCTAAATAGATACTCGCAAAGGCGACTGATATATCCTTAGGTCGAAGAGATACCAGAGAAATAGTAGGTCTTAGGCTAGTCCTTGTCCAGTA$TACTGGACAAGGACTAGCCTAAGACCTACTATTTCTCTGGTATCTCTTCGACCTAAGGATATATCAGTCGCCTTTGCGAGTATCTATTTAGGTGTCCCTA$ACGCACCCCGGCATTCGTCGACTCTACACTTAGTGGAACATACAAATTCGCTCGCAGGAGCGCCTCATACATTCTAACGCAGTGATCTTCGGCTGAGACT$AGTCTCAGCCGAAGATCACTGCGTTAGAATGTATGAGGCGCTCCTGCGAGCGAATTTGTATGTTCCACTAAGTGTAGAGTCGACGAATGCCGGGGTGCGT$"[..], "complex sentinels"),
             ];
        let num_rand = 100;
        let rand_cases: Vec<_> = (0..num_rand).map(|i| rand_seqs(10, i * 10)).collect();
        test_cases.extend(
            rand_cases
                .iter()
                .map(|case| (case.as_ref(), "rand test case")),
        );

        for &(text, test_name) in test_cases.iter() {
            for &sample_rate in &[2, 3, 5, 16] {
                let alphabet = dna::n_alphabet();
                let sa = suffix_array(text);
                let bwt = bwt(text, &sa);
                let less = less(&bwt, &alphabet);
                let occ = Occ::new(&bwt, 3, &alphabet);
                let sampled = sa.sample(text, &bwt, &less, &occ, sample_rate);

                for i in 0..sa.len() {
                    let sa_idx = sa.get(i).unwrap();
                    let sampled_idx = sampled.get(i).unwrap();
                    assert_eq!(
                        sa_idx,
                        sampled_idx,
                        "Failed:\n{}\n{}\nat index {} do not match in test: {} (sample rate: {})",
                        str::from_utf8(&text[sa_idx..]).unwrap(),
                        str::from_utf8(&text[sampled_idx..]).unwrap(),
                        i,
                        test_name,
                        sample_rate
                    );
                }
            }
        }
    }

    #[test]
    fn can_construct_sa_for_usize() {
        let text: Vec<usize> = vec![3, 2, 2, 4, 4, 1, 2, 1, 0];
        let sa = suffix_array_int(&text);
        assert_eq!(sa, vec![8, 7, 5, 6, 1, 2, 0, 4, 3]);
    }
}

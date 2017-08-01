// Copyright 2014-2015 Johannes Köster, Vadim Nazarov, Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Banded Smith-Waterman alignment for fast comparison of long strings.
//! Use sparse dynamic programming to find a 'backbone' alignment from exact
//! k-mer matches, then compute the SW alignment in a 'band' surrounding the
//! backbone, with a configurable width w. This method is not guaranteed
//! to recover the Smith-Waterman alignment, but will usually find the same
//! alignment if a) there is a reasonable density of exact k-mer matches
//! between the sequences, and b) the width parameter w is larger than the
//! excursion of the alignment path from diagonal between successive kmer
//! matches.  This technique is employed in long-read aligners (e.g. BLASR and BWA)
//! to drastically reduce runtime compared to Smith Waterman.
//! Complexity roughly O(min(m,n) * w)
//!
//! # Example
//!
//! ```
//! use bio::alignment::pairwise::banded::*;
//! use bio::alignment::sparse::hash_kmers;
//! use bio::alignment::pairwise::{MIN_SCORE, Scoring};
//! use bio::alignment::AlignmentOperation::*;
//! use std::iter::repeat;
//!
//! let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
//! let y = b"AGCACACGTGTGCGCTATACAGTACACGTGTCACAGTTGTACTAGCATGAC";
//! let score = |a: u8, b: u8| if a == b {1i32} else {-1i32};
//! let k = 8;  // kmer match length
//! let w = 6;  // Window size for creating the band
//! let mut aligner = Aligner::new(-5, -1, score, k, w);
//! let alignment = aligner.local(x, y);
//! // aligner.global(x, y), aligner.semiglobal(x, y) are also supported
//! assert_eq!(alignment.ystart, 0);
//! assert_eq!(alignment.xstart, 0);
//!
//! // For cases where the reference is reused multiple times, we can invoke the
//! // pre-hashed version of the solver
//! let x = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCA";
//! let y = b"CAGTTGTACTAGCATGACCAGTTGTACTAGCATGACAGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCA\
//!     CAGTTGTACTAGCATGACCAGTTGTACTAGCATGAC";
//! let y_kmers_hash = hash_kmers(y, k);
//! let alignment = aligner.semiglobal_with_prehash(x, y, &y_kmers_hash);
//! assert_eq!(alignment.score, 37);
//!
//! // In addition to the standard modes (Global, Semiglobal and Local), a custom alignment
//! // mode is supported which supports a user-specified clipping penalty. Clipping is a
//! // special boundary condition where you are allowed to clip off the beginning/end of
//! // the sequence for a fixed penalty. See bio::alignment::pairwise for a more detailed
//! // explanation
//!
//! // The following example considers a modification of the semiglobal mode where you are allowed
//! // to skip a prefix of the target sequence x, for a penalty of -10, but you have to consume
//! // the rest of the string in the alignment
//!
//! let scoring = Scoring {
//!     gap_open: -5,
//!     gap_extend: -1,
//!     match_fn: |a: u8, b: u8| if a == b {1i32} else {-3i32},
//!     xclip_prefix: -10,
//!     xclip_suffix: MIN_SCORE,
//!     yclip_prefix: 0,
//!     yclip_suffix: 0
//! };
//! let x = b"GGGGGGACGTACGTACGTGTGCATCATCATGTGCGTATCATAGATAGATGTAGATGATCCACAGT";
//! let y = b"AAAAACGTACGTACGTGTGCATCATCATGTGCGTATCATAGATAGATGTAGATGATCCACAGTAAAA";
//! let mut aligner = Aligner::with_capacity_and_scoring(x.len(), y.len(), scoring, k, w);
//! let alignment = aligner.custom(x, y);
//! println!("{}", alignment.pretty(x,y));
//! assert_eq!(alignment.score, 49);
//! let mut correct_ops = Vec::new();
//! correct_ops.push(Yclip(4));
//! correct_ops.push(Xclip(6));
//! correct_ops.extend(repeat(Match).take(59));
//! correct_ops.push(Yclip(4));
//! assert_eq!(alignment.operations, correct_ops);
//!
//! // aligner.custom_with_prehash(x, y, &y_kmers_hash) is also supported
//! ```

use std::i32;
use alignment::{Alignment, AlignmentOperation};
use utils::TextSlice;
use std::cmp::min;
use std::cmp::max;
use std::ops::Range;

use super::*;
use alignment::sparse;
use alignment::sparse::HashMapFx;
use alignment::pairwise::Scoring;

const MAX_CELLS: usize = 100000;

/// A banded implementation of Smith-Waterman aligner (SWA).
/// Unlike the full SWA, this implementation computes the alignment between a pair of sequences
/// only inside a 'band' withing the dynamic programming matrix. The band is constructed using the
/// Sparse DP routie (see sparse::sdpkpp), which uses kmer matches to build the best common
/// subsequence (including gap penalties) between the two strings. The band is constructed around
/// this subsequence (using the window length 'w'), filling in the gaps.
///
/// In the case where there are no k-mer matches, the  aligner will fall back to a full alignment,
/// by setting the band to contain the full matrix.
///
/// Banded aligner will proceed to compute the alignment only when the total number of cells
/// in the band is less than MAX_CELLS (currently set to 10 million), otherwise it returns an
/// empty alignment
#[allow(non_snake_case)]
pub struct Aligner<F: MatchFunc> {
    S: [Vec<i32>; 2],
    I: [Vec<i32>; 2],
    D: [Vec<i32>; 2],
    Lx: Vec<usize>,
    Ly: Vec<usize>,
    Sn: Vec<i32>,
    traceback: Traceback,
    scoring: Scoring<F>,

    band: Band,
    k: usize,
    w: usize,
}


const DEFAULT_ALIGNER_CAPACITY: usize = 200;

impl<F: MatchFunc> Aligner<F> {
    /// Create new aligner instance with given gap open and gap extend penalties
    /// and the score function.
    ///
    /// # Arguments
    ///
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    /// * `k` - kmer length used in constructing the band
    /// * `w` - width of the band
    ///
    pub fn new(gap_open: i32, gap_extend: i32, match_fn: F, k: usize, w: usize) -> Self {
        Aligner::with_capacity(DEFAULT_ALIGNER_CAPACITY,
                               DEFAULT_ALIGNER_CAPACITY,
                               gap_open,
                               gap_extend,
                               match_fn,
                               k,
                               w)
    }

    /// Create new aligner instance. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `gap_open` - the score for opening a gap (should be negative)
    /// * `gap_extend` - the score for extending a gap (should be negative)
    /// * `match_fn` - function that returns the score for substitutions (also see bio::scores)
    /// * `k` - kmer length used in constructing the band
    /// * `w` - width of the band
    ///
    pub fn with_capacity(m: usize,
                         n: usize,
                         gap_open: i32,
                         gap_extend: i32,
                         match_fn: F,
                         k: usize,
                         w: usize)
                         -> Self {

        Aligner {
            band: Band::new(m, n),
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring: Scoring::new(gap_open, gap_extend, match_fn),
            k: k,
            w: w,
        }
    }

    /// Create new aligner instance with scoring and size hint. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `scoring` - the scoring struct
    /// * `k` - kmer length used in constructing the band
    /// * `w` - width of the band
    ///
    pub fn with_capacity_and_scoring(m: usize,
                                     n: usize,
                                     scoring: Scoring<F>,
                                     k: usize,
                                     w: usize)
                                     -> Self {

        assert!(scoring.gap_open <= 0, "gap_open can't be positive");
        assert!(scoring.gap_extend <= 0, "gap_extend can't be positive");
        assert!(scoring.xclip_prefix <= 0,
                "Clipping penalty (x prefix) can't be positive");
        assert!(scoring.xclip_suffix <= 0,
                "Clipping penalty (x suffix) can't be positive");
        assert!(scoring.yclip_prefix <= 0,
                "Clipping penalty (y prefix) can't be positive");
        assert!(scoring.yclip_suffix <= 0,
                "Clipping penalty (y suffix) can't be positive");

        Aligner {
            band: Band::new(m, n),
            S: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            I: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            D: [Vec::with_capacity(m + 1), Vec::with_capacity(m + 1)],
            Lx: Vec::with_capacity(n + 1),
            Ly: Vec::with_capacity(m + 1),
            Sn: Vec::with_capacity(m + 1),
            traceback: Traceback::with_capacity(m, n),
            scoring: scoring,
            k: k,
            w: w,
        }
    }

    /// Create new aligner instance with scoring and size hint. The size hints help to
    /// avoid unnecessary memory allocations.
    ///
    /// # Arguments
    ///
    /// * `m` - the expected size of x
    /// * `n` - the expected size of y
    /// * `scoring` - the scoring struct
    /// * `k` - kmer length used in constructing the band
    /// * `w` - width of the band
    ///
    pub fn with_scoring(scoring: Scoring<F>, k: usize, w: usize) -> Self {

        Aligner::with_capacity_and_scoring(DEFAULT_ALIGNER_CAPACITY,
                                           DEFAULT_ALIGNER_CAPACITY,
                                           scoring,
                                           k,
                                           w)
    }

    /// Compute the alignment with custom clip penalties
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    ///
    pub fn custom(&mut self, x: TextSlice, y: TextSlice) -> Alignment {
        self.band = Band::create(x, y, self.k, self.w, &self.scoring);
        self.compute_alignment(x, y)
    }

    /// Compute the alignment with custom clip penalties with 'y' being pre-hashed
    /// (see sparse::hash_kmers)
    ///
    /// # Arguments
    ///
    /// * `x` - Textslice
    /// * `y` - Textslice
    ///
    pub fn custom_with_prehash(&mut self,
                               x: TextSlice,
                               y: TextSlice,
                               y_kmer_hash: &HashMapFx<&[u8], Vec<u32>>)
                               -> Alignment {
        self.band = Band::create_with_prehash(x, y, self.k, self.w, &self.scoring, y_kmer_hash);
        self.compute_alignment(x, y)
    }

    // Computes the alignment. The band needs to be populated prior
    // to calling this function
    fn compute_alignment(&mut self, x: TextSlice, y: TextSlice) -> Alignment {

        if self.band.num_cells() > MAX_CELLS {
            // Too many cells in the band. Return an empty alignment
            return Alignment {
                       score: MIN_SCORE,
                       ystart: 0,
                       xstart: 0,
                       yend: 0,
                       xend: 0,
                       ylen: 0,
                       xlen: 0,
                       operations: Vec::new(),
                       mode: AlignmentMode::Custom,
                   };
        }

        let (m, n) = (x.len(), y.len());
        self.traceback.init(m, n);

        for k in 0..2 {
            self.I[k].clear();
            self.D[k].clear();
            self.S[k].clear();
            self.D[k].extend(repeat(MIN_SCORE).take(m + 1));
            self.I[k].extend(repeat(MIN_SCORE).take(m + 1));
            self.S[k].extend(repeat(MIN_SCORE).take(m + 1));
        }
        self.Lx.clear();
        self.Lx.extend(repeat(0usize).take(n + 1));
        self.Ly.clear();
        self.Ly.extend(repeat(0usize).take(m + 1));
        self.Sn.clear();
        self.Sn.extend(repeat(MIN_SCORE).take(m + 1));

        {
            // Handle j = 0
            let curr = 0;
            let i_start = self.band.ranges[0].start;
            let i_end = self.band.ranges[0].end;
            if i_start == 0 {
                self.S[curr][0] = 0;
            }

            for i in max(1, i_start)..i_end {
                let mut tb = TracebackCell::new();
                tb.set_all(TB_START);
                if i == 1 {
                    self.I[curr][i] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_i_bits(TB_START);
                } else {
                    // Insert all i characters
                    let i_score = self.scoring.gap_open + self.scoring.gap_extend * (i as i32);
                    let c_score = self.scoring.xclip_prefix + self.scoring.gap_open +
                                  self.scoring.gap_extend; // Clip then insert
                    if i_score > c_score {
                        self.I[curr][i] = i_score;
                        tb.set_i_bits(TB_INS);
                    } else {
                        self.I[curr][i] = c_score;
                        tb.set_i_bits(TB_XCLIP_PREFIX);
                    }
                }

                if i == m {
                    tb.set_s_bits(TB_XCLIP_SUFFIX);
                }

                if self.I[curr][i] > self.S[curr][i] {
                    self.S[curr][i] = self.I[curr][i];
                    tb.set_s_bits(TB_INS);
                }

                if self.scoring.xclip_prefix > self.S[curr][i] {
                    self.S[curr][i] = self.scoring.xclip_prefix;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                // Track the score if we do a suffix clip (x) after this character
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[0] = m - i;
                }

                self.traceback.set(i, 0, tb);
            }

            for i in i_end..min(m + 1, self.band.ranges[min(n, 1)].end) {
                self.S[curr][i] = MIN_SCORE;
                self.I[curr][i] = MIN_SCORE;
            }

            if i_end < (m + 1) {
                self.S[curr][m] = MIN_SCORE;
            }
        }

        for j in 1..n + 1 {
            let curr = j % 2;
            let prev = 1 - curr;

            let i_start = self.band.ranges[j].start;
            let i_end = self.band.ranges[j].end;

            if i_start == 0 {
                // Handle i = 0
                let mut tb = TracebackCell::new();
                self.I[curr][0] = MIN_SCORE;

                if j == 1 {
                    self.D[curr][0] = self.scoring.gap_open + self.scoring.gap_extend;
                    tb.set_d_bits(TB_START);
                } else {
                    // Delete all j characters
                    let d_score = self.scoring.gap_open + self.scoring.gap_extend * (j as i32);
                    let c_score = self.scoring.yclip_prefix + self.scoring.gap_open +
                                  self.scoring.gap_extend;
                    if d_score > c_score {
                        self.D[curr][0] = d_score;
                        tb.set_d_bits(TB_DEL);
                    } else {
                        self.D[curr][0] = c_score;
                        tb.set_d_bits(TB_YCLIP_PREFIX);
                    }
                }

                if self.D[curr][0] > self.scoring.yclip_prefix {
                    self.S[curr][0] = self.D[curr][0];
                    tb.set_s_bits(TB_DEL);
                } else {
                    self.S[curr][0] = self.scoring.yclip_prefix;
                    tb.set_s_bits(TB_YCLIP_PREFIX);
                }

                // Track the score if we do suffix clip (y) from here
                if self.S[curr][0] + self.scoring.yclip_suffix > self.Sn[0] {
                    self.Sn[0] = self.S[curr][0] + self.scoring.yclip_suffix;
                    self.Ly[0] = n - j;
                }
                self.traceback.set(0, j, tb);
            }

            for i in i_start.saturating_sub(1)..i_start {
                self.S[curr][i] = MIN_SCORE;
                self.I[curr][i] = MIN_SCORE;
                self.D[curr][i] = MIN_SCORE;
            }
            self.S[curr][m] = MIN_SCORE;

            let q = y[j - 1];
            let xclip_score = self.scoring.xclip_prefix +
                              max(self.scoring.yclip_prefix,
                                  self.scoring.gap_open + self.scoring.gap_extend * (j as i32));

            for i in max(1, i_start)..i_end {

                let p = x[i - 1];
                let mut tb = TracebackCell::new();

                let m_score = self.S[prev][i - 1] + self.scoring.match_fn.score(p, q);

                let i_score = self.I[curr][i - 1] + self.scoring.gap_extend;
                let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
                let best_i_score;
                if i_score > s_score {
                    best_i_score = i_score;
                    tb.set_i_bits(TB_INS);
                } else {
                    best_i_score = s_score;
                    tb.set_i_bits(self.traceback.get(i - 1, j).get_s_bits());
                }

                let d_score = self.D[prev][i] + self.scoring.gap_extend;
                let s_score = self.S[prev][i] + self.scoring.gap_open + self.scoring.gap_extend;
                let best_d_score;
                if d_score > s_score {
                    best_d_score = d_score;
                    tb.set_d_bits(TB_DEL);
                } else {
                    best_d_score = s_score;
                    tb.set_d_bits(self.traceback.get(i, j - 1).get_s_bits());
                }

                if i == m {
                    tb.set_s_bits(TB_XCLIP_SUFFIX);
                } else {
                    self.S[curr][i] = MIN_SCORE;
                }
                let mut best_s_score = self.S[curr][i];

                if m_score > best_s_score {
                    best_s_score = m_score;
                    tb.set_s_bits(if p == q { TB_MATCH } else { TB_SUBST });
                }

                if best_i_score > best_s_score {
                    best_s_score = best_i_score;
                    tb.set_s_bits(TB_INS);
                }

                if best_d_score > best_s_score {
                    best_s_score = best_d_score;
                    tb.set_s_bits(TB_DEL);
                }

                if xclip_score > best_s_score {
                    best_s_score = xclip_score;
                    tb.set_s_bits(TB_XCLIP_PREFIX);
                }

                let yclip_score = self.scoring.yclip_prefix + self.scoring.gap_open +
                                  self.scoring.gap_extend * (i as i32);
                if yclip_score > best_s_score {
                    best_s_score = yclip_score;
                    tb.set_s_bits(TB_YCLIP_PREFIX);
                }

                self.S[curr][i] = best_s_score;
                self.I[curr][i] = best_i_score;
                self.D[curr][i] = best_d_score;

                // Track the score if we do suffix clip (x) from here
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[j] = m - i;
                }

                // Track the score if we do suffix clip (y) from here
                if self.S[curr][i] + self.scoring.yclip_suffix > self.Sn[i] {
                    self.Sn[i] = self.S[curr][i] + self.scoring.yclip_suffix;
                    self.Ly[i] = n - j;
                }

                self.traceback.set(i, j, tb);

            }

            // Suffix clip (y) from i = m and reset S[curr][m] if required
            if self.S[curr][m] + self.scoring.yclip_suffix > self.Sn[m] {
                self.Sn[m] = self.S[curr][m] + self.scoring.yclip_suffix;
                self.Ly[m] = n - j;
            }
            if i_end < (m + 1) {
                self.traceback.get_mut(m, j).set_s_bits(TB_XCLIP_SUFFIX);
                self.S[curr][m] = MIN_SCORE;
            }

            for i in i_end..min(m + 1, self.band.ranges[min(n, j + 1)].end) {
                self.S[curr][i] = MIN_SCORE;
                self.I[curr][i] = MIN_SCORE;
                self.D[curr][i] = MIN_SCORE;
            }

        }

        // Handle suffix clipping in the j=n case
        for i in 0..m + 1 {
            let j = n;
            let curr = j % 2;
            if self.Sn[i] > self.S[curr][i] {
                self.S[curr][i] = self.Sn[i];
                self.traceback.get_mut(i, j).set_s_bits(TB_YCLIP_SUFFIX);
            }
            if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                self.Lx[j] = m - i;
                self.traceback.get_mut(m, j).set_s_bits(TB_XCLIP_SUFFIX);
            }
        }

        // Since there could be a change in the last column of S,
        // recompute the last colum of I as this could also change
        for i in max(1, self.band.ranges[n].start)..self.band.ranges[n].end {
            let j = n;
            let curr = j % 2;
            let s_score = self.S[curr][i - 1] + self.scoring.gap_open + self.scoring.gap_extend;
            if s_score > self.I[curr][i] {
                self.I[curr][i] = s_score;
                let s_bit = self.traceback.get(i - 1, j).get_s_bits();
                self.traceback.get_mut(i, j).set_i_bits(s_bit);
            }
            if s_score > self.S[curr][i] {
                self.S[curr][i] = s_score;
                self.traceback.get_mut(i, j).set_s_bits(TB_INS);
                if self.S[curr][i] + self.scoring.xclip_suffix > self.S[curr][m] {
                    self.S[curr][m] = self.S[curr][i] + self.scoring.xclip_suffix;
                    self.Lx[j] = m - i;
                    self.traceback.get_mut(m, j).set_s_bits(TB_XCLIP_SUFFIX);
                }
            }
        }

        let mut i = m;
        let mut j = n;
        let mut ops = Vec::with_capacity(x.len());
        let mut xstart: usize = 0usize;
        let mut ystart: usize = 0usize;
        let mut xend = m;
        let mut yend = n;

        let mut last_layer = self.traceback.get(i, j).get_s_bits();

        loop {
            let next_layer: u16;
            match last_layer {
                TB_START => break,
                TB_INS => {
                    ops.push(AlignmentOperation::Ins);
                    next_layer = self.traceback.get(i, j).get_i_bits();
                    i -= 1;
                }
                TB_DEL => {
                    ops.push(AlignmentOperation::Del);
                    next_layer = self.traceback.get(i, j).get_d_bits();
                    j -= 1;
                }
                TB_MATCH => {
                    ops.push(AlignmentOperation::Match);
                    next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_SUBST => {
                    ops.push(AlignmentOperation::Subst);
                    next_layer = self.traceback.get(i - 1, j - 1).get_s_bits();
                    i -= 1;
                    j -= 1;
                }
                TB_XCLIP_PREFIX => {
                    ops.push(AlignmentOperation::Xclip(i));
                    xstart = i;
                    i = 0;
                    next_layer = self.traceback.get(0, j).get_s_bits();
                }
                TB_XCLIP_SUFFIX => {
                    ops.push(AlignmentOperation::Xclip(self.Lx[j]));
                    i -= self.Lx[j];
                    xend = i;
                    next_layer = self.traceback.get(i, j).get_s_bits();
                }
                TB_YCLIP_PREFIX => {
                    ops.push(AlignmentOperation::Yclip(j));
                    ystart = j;
                    j = 0;
                    next_layer = self.traceback.get(i, 0).get_s_bits();
                }
                TB_YCLIP_SUFFIX => {
                    ops.push(AlignmentOperation::Yclip(self.Ly[i]));
                    j -= self.Ly[i];
                    yend = j;
                    next_layer = self.traceback.get(i, j).get_s_bits();
                }
                _ => panic!("Dint expect this!"),
            }
            last_layer = next_layer;
            // println!("{} of {}, {} of {} - {}", i, m, j, n, last_layer);
        }

        // Handle the case when the traceback ends outside the band other than at (0, 0)
        if i != 0 {
            ops.push(AlignmentOperation::Xclip(i));
            xstart = i;
        }
        if j != 0 {
            ops.push(AlignmentOperation::Yclip(j));
            ystart = j;
        }

        ops.reverse();
        Alignment {
            score: self.S[n % 2][m],
            ystart: ystart,
            xstart: xstart,
            yend: yend,
            xend: xend,
            ylen: n,
            xlen: m,
            operations: ops,
            mode: AlignmentMode::Custom,
        }
    }

    /// Calculate global alignment of x against y.
    pub fn global(&mut self, x: TextSlice, y: TextSlice) -> Alignment {

        // Store the current clip penalties
        let clip_penalties = [self.scoring.xclip_prefix,
                              self.scoring.xclip_suffix,
                              self.scoring.yclip_prefix,
                              self.scoring.yclip_suffix];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = MIN_SCORE;
        self.scoring.yclip_suffix = MIN_SCORE;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Global;

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    pub fn semiglobal(&mut self, x: TextSlice, y: TextSlice) -> Alignment {

        // Store the current clip penalties
        let clip_penalties = [self.scoring.xclip_prefix,
                              self.scoring.xclip_suffix,
                              self.scoring.yclip_prefix,
                              self.scoring.yclip_suffix];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Xclip and Yclip from alignment.operations
        alignment.filter_clip_operations();

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate semiglobal alignment of x against y (x is global, y is local).
    /// This function accepts the hash map of the kmers of y. This is useful
    /// in cases where we are interested in repeated alignment of different
    /// queries against the same reference. The user can precompute the HashMap
    /// using sparse::hash_kmers and invoke this function to speed up the
    /// alignment computation.
    pub fn semiglobal_with_prehash(&mut self,
                                   x: TextSlice,
                                   y: TextSlice,
                                   y_kmer_hash: &HashMapFx<&[u8], Vec<u32>>)
                                   -> Alignment {

        // Store the current clip penalties
        let clip_penalties = [self.scoring.xclip_prefix,
                              self.scoring.xclip_suffix,
                              self.scoring.yclip_prefix,
                              self.scoring.yclip_suffix];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = MIN_SCORE;
        self.scoring.xclip_suffix = MIN_SCORE;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom_with_prehash(x, y, y_kmer_hash);
        alignment.mode = AlignmentMode::Semiglobal;

        // Filter out Xclip and Yclip from alignment.operations
        alignment.filter_clip_operations();

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    /// Calculate local alignment of x against y.
    pub fn local(&mut self, x: TextSlice, y: TextSlice) -> Alignment {

        // Store the current clip penalties
        let clip_penalties = [self.scoring.xclip_prefix,
                              self.scoring.xclip_suffix,
                              self.scoring.yclip_prefix,
                              self.scoring.yclip_suffix];

        // Temporarily Over-write the clip penalties
        self.scoring.xclip_prefix = 0;
        self.scoring.xclip_suffix = 0;
        self.scoring.yclip_prefix = 0;
        self.scoring.yclip_suffix = 0;

        // Compute the alignment
        let mut alignment = self.custom(x, y);
        alignment.mode = AlignmentMode::Local;

        // Filter out Xclip and Yclip from alignment.operations
        alignment.filter_clip_operations();

        // Set the clip penalties to the original values
        self.scoring.xclip_prefix = clip_penalties[0];
        self.scoring.xclip_suffix = clip_penalties[1];
        self.scoring.yclip_prefix = clip_penalties[2];
        self.scoring.yclip_suffix = clip_penalties[3];

        alignment
    }

    #[allow(dead_code)]
    fn visualize(&self, alignment: &Alignment) {
        // First populate the band
        let mut view = vec!['.'; self.band.rows * self.band.cols];
        let index = |i, j| i * self.band.cols + j;
        for j in 0..self.band.ranges.len() {
            let range = &self.band.ranges[j];
            for i in range.start..range.end {
                view[index(i, j)] = 'x';
            }
        }

        // Populate the path
        let path = alignment.path();
        for p in path {
            view[index(p.0, p.1)] = '\\';
        }

        for i in 0..self.band.rows {
            for j in 0..self.band.cols {
                print!("{}", view[index(i, j)]);
            }
            println!("");
        }
    }
}

trait MatchPair {
    fn continues(&self, p: Option<(u32, u32)>) -> bool;
}

impl MatchPair for (u32, u32) {
    fn continues(&self, p: Option<(u32, u32)>) -> bool {
        match p {
            Some(_p) => (self.0 == _p.0 + 1 && self.1 == _p.1 + 1),
            None => false,
        }
    }
}

#[derive(Clone, Debug)]
struct Band {
    rows: usize,
    cols: usize,
    ranges: Vec<Range<usize>>,
}

impl Band {
    // Create new Band instance with given size
    //
    // # Arguments
    //
    // * `m` - the expected size of x
    // * `n` - the expected size of y
    //
    fn new(m: usize, n: usize) -> Self {

        let mut ranges: Vec<Range<usize>> = Vec::with_capacity(n + 1);
        for _ in 0..n + 1 {
            ranges.push(m + 1..0);
        }

        Band {
            rows: m + 1,
            cols: n + 1,
            ranges: ranges,
        }

    }

    // Add cells around a kmer of length 'k', starting at 'start', which are within a
    // distance of 'w' in x or y directions to the band.
    fn add_kmer(&mut self, start: (u32, u32), k: usize, w: usize) {

        let (r, c) = (start.0 as usize, start.1 as usize);
        // println!("{} {} {}", r, k, self.rows);
        debug_assert!(r + k <= self.rows);
        debug_assert!(c + k <= self.cols);

        if k == 0 {
            return;
        }

        let i = r.saturating_sub(w);
        for j in c.saturating_sub(w)..min(c + w + 1, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, i);
        }

        let mut i = r.saturating_sub(w);
        for j in min(c + w, self.cols)..min(c + k + w, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, i);
            i += 1;
        }

        let mut i = r + w + k;
        let mut j = (c + k - 1).saturating_sub(w);
        loop {
            if j <= c.saturating_sub(w) {
                break;
            }
            j -= 1;
            i -= 1;
            self.ranges[j].end = max(self.ranges[j].end, min(i, self.rows));
        }

        let i = min(r + w + k, self.rows);
        for j in (c + k - 1).saturating_sub(w)..min(c + k + w, self.cols) {
            self.ranges[j].end = max(self.ranges[j].end, i);
        }

    }

    // Add cells around a specific position to the band. An cell which is within 'w' distance
    // in x or y directions are added
    fn add_entry(&mut self, pos: (u32, u32), w: usize) {
        let (r, c) = (pos.0 as usize, pos.1 as usize);

        let istart = r.saturating_sub(w);
        let iend = min(r + w + 1, self.rows);
        for j in c.saturating_sub(w)..min(c + w + 1, self.cols) {
            self.ranges[j].start = min(self.ranges[j].start, istart);
            self.ranges[j].end = max(self.ranges[j].end, iend);
        }
    }

    // Each gap generates a line from the start to end.
    fn add_gap(&mut self, start: (u32, u32), end: (u32, u32), w: usize) {
        let nrows = end.0 - start.0;
        let ncols = end.1 - start.1;
        if nrows > ncols {
            for r in start.0..end.0 {
                let c = start.1 + (end.1 - start.1) * (r - start.0) / (end.0 - start.0);
                self.add_entry((r, c), w);
            }
        } else {
            for c in start.1..end.1 {
                let r = start.0 + (end.0 - start.0) * (c - start.1) / (end.1 - start.1);
                self.add_entry((r, c), w);
            }
        }
    }


    // The band needs to start either at (0,0) or at a point that is zero score from (0,0).
    // This naturally sets the start positions correctly for global, semiglobal and local
    // modes. Similarly the band has to either end at (m,n) or at a point from which there is
    // a zero score path to (m,n).
    //
    // At the minimum, irrespective of the score (0,0)->start or end->(m,n), we extend the band
    // diagonally for a length "lazy_extend"(2k) or when it hits the corner, whichever happens first
    //
    // start - the index of the first matching kmer in LCSk++
    // end - the index of the last matching kmer in LCSk++
    //
    fn set_boundaries<F: MatchFunc>(&mut self,
                                    start: (u32, u32),
                                    end: (u32, u32),
                                    k: usize,
                                    w: usize,
                                    scoring: &Scoring<F>) {

        let lazy_extend: usize = 2 * k;

        // -------------- START --------------
        // Nothing to do if the start is already at (0,0)
        let (r, c) = (start.0 as usize, start.1 as usize);
        if !(r == 0usize && c == 0usize) {
            let mut score_to_start = if r > 0 { scoring.xclip_prefix } else { 0i32 };
            score_to_start += if c > 0 { scoring.yclip_prefix } else { 0i32 };

            if score_to_start == 0 {
                // Just do a "lazy_extend"
                // First diagonally
                let d = min(lazy_extend, min(r, c));
                self.add_kmer(((r - d) as u32, (c - d) as u32), d, w);

                // If we hit one of the edges before completing lazy_extend
                self.add_gap((r.saturating_sub(lazy_extend) as u32,
                              c.saturating_sub(lazy_extend) as u32),
                             ((r - d) as u32, (c - d) as u32),
                             w);

            } else {
                // we need to find a zero cost cell

                // First try the diagonal
                let diagonal_score;
                if r > c {
                    // We will hit (r-c, 0)
                    diagonal_score = scoring.xclip_prefix;
                } else if c > r {
                    // We will hit (0, c-r)
                    diagonal_score = scoring.yclip_prefix;
                } else {
                    diagonal_score = 0;
                }

                if diagonal_score == 0 {
                    let d = min(r, c);
                    self.add_kmer(((r - d) as u32, (c - d) as u32), d, w);
                    // Make sure we do at least "lazy_extend" extension
                    let start = (r.saturating_sub(lazy_extend) as u32,
                                 c.saturating_sub(lazy_extend) as u32);
                    let end = ((r - d) as u32, (c - d) as u32);
                    if (start.0 <= end.0) && (start.1 <= end.1) {
                        self.add_gap(start, end, w);
                    }
                } else {
                    // Band to origin
                    self.add_gap((0u32, 0u32), start, w);
                }
            }
        }


        // -------------- END --------------
        // Nothing to do if the last kmer ends at (m, n)
        let (r, c) = (end.0 as usize + k, end.1 as usize + k);
        debug_assert!(r <= self.rows);
        debug_assert!(c <= self.cols);
        if !(r == self.rows && c == self.cols) {
            let mut score_from_end = if r == self.rows {
                0
            } else {
                scoring.xclip_suffix
            };
            score_from_end += if c == self.cols {
                0
            } else {
                scoring.yclip_suffix
            };

            if score_from_end == 0 {
                // Just a lazy_extend
                let d = min(lazy_extend, min(self.rows - r, self.cols - c));
                self.add_kmer((r as u32, c as u32), d, w);

                let r1 = min(self.rows, r + d) - 1;
                let c1 = min(self.cols, c + d) - 1;
                let r2 = min(self.rows, r + lazy_extend);
                let c2 = min(self.cols, c + lazy_extend);
                if (r1 <= r2) && (c1 <= c2) {
                    self.add_gap((r1 as u32, c1 as u32), (r2 as u32, c2 as u32), w);
                }

            } else {
                // we need to find a zero cost cell

                // First try the diagonal
                let dr = self.rows - r;
                let dc = self.cols - c;
                let diagonal_score;
                if dr > dc {
                    // We will hit (r+dc, self.cols)
                    diagonal_score = scoring.xclip_suffix;
                } else if dc > dr {
                    // We will hit (self.rows, c+dr)
                    diagonal_score = scoring.yclip_suffix;
                } else {
                    // We will hit the corner
                    diagonal_score = 0;
                }

                if diagonal_score == 0 {
                    let d = min(dr, dc);
                    self.add_kmer((r as u32, c as u32), d, w);
                    // Make sure we do at least "lazy_extend" extension
                    let r1 = min(self.rows, r + d) - 1;
                    let c1 = min(self.cols, c + d) - 1;
                    let r2 = min(self.rows, r + lazy_extend);
                    let c2 = min(self.cols, c + lazy_extend);
                    if (r1 <= r2) && (c1 <= c2) {
                        self.add_gap((r1 as u32, c1 as u32), (r2 as u32, c2 as u32), w);
                    }
                } else {
                    // Band to lower right corner
                    let rows = self.rows as u32;
                    let cols = self.cols as u32;
                    self.add_gap((r as u32, c as u32), (rows as u32, cols as u32), w);
                }

            }
        }
    }


    fn create<F: MatchFunc>(x: TextSlice,
                            y: TextSlice,
                            k: usize,
                            w: usize,
                            scoring: &Scoring<F>)
                            -> Band {

        let matches = sparse::find_kmer_matches(x, y, k);
        Band::create_with_matches(x, y, k, w, scoring, matches)
    }

    fn create_with_prehash<F: MatchFunc>(x: TextSlice,
                                         y: TextSlice,
                                         k: usize,
                                         w: usize,
                                         scoring: &Scoring<F>,
                                         y_kmer_hash: &HashMapFx<&[u8], Vec<u32>>)
                                         -> Band {

        let matches = sparse::find_kmer_matches_seq2_hashed(x, y_kmer_hash, k);
        Band::create_with_matches(x, y, k, w, scoring, matches)
    }

    fn create_with_matches<F: MatchFunc>(x: TextSlice,
                                         y: TextSlice,
                                         k: usize,
                                         w: usize,
                                         scoring: &Scoring<F>,
                                         matches: Vec<(u32, u32)>)
                                         -> Band {

        let mut band = Band::new(x.len(), y.len());

        if matches.len() == 0 {
            band.full_matrix();
            return band;
        }
        let res = sparse::sdpkpp(&matches, k, 2, scoring.gap_open, scoring.gap_extend);
        let ps = res.path[0];
        let pe = res.path[res.path.len() - 1];

        // Set the boundaries
        band.set_boundaries(matches[ps], matches[pe], k, w, scoring);

        // for idx in &res.path {
        //     println!("{:?}", matches[*idx]);
        // }

        // println!("sparse: rstart:{} tstart:{} rend:{}, tend:{}, hits:{}",
        //          matches[ps].0,
        //          matches[ps].1,
        //          matches[pe].0,
        //          matches[pe].1,
        //          res.score);

        let mut prev: Option<(u32, u32)> = None;

        for idx in res.path {
            let curr = matches[idx];
            if curr.continues(prev) {
                let p = prev.unwrap();
                band.add_entry((p.0 + k as u32, p.1 + k as u32), w);
            } else {
                match prev {
                    Some(p) => band.add_gap((p.0 + (k - 1) as u32, p.1 + (k - 1) as u32), curr, w),
                    _ => {}
                }
                band.add_kmer(curr, k, w);
            }
            prev = Some(curr);
        }
        band
    }

    fn full_matrix(&mut self) {
        self.ranges.clear();
        for _ in 0..self.cols {
            self.ranges.push(0..self.rows);
        }
    }

    fn num_cells(&self) -> usize {
        let mut banded_cells = 0;
        for j in 0..self.ranges.len() {
            banded_cells += self.ranges[j].end.saturating_sub(self.ranges[j].start);
        }
        banded_cells
    }

    #[allow(dead_code)]
    fn visualize(&self) {
        let mut view = vec!['.'; self.rows * self.cols];
        let index = |i, j| i * self.cols + j;
        for j in 0..self.ranges.len() {
            let range = &self.ranges[j];
            for i in range.start..range.end {
                view[index(i, j)] = 'x';
            }
        }

        for i in 0..self.rows {
            for j in 0..self.cols {
                print!("{}", view[index(i, j)]);
            }
            println!("");
        }
    }

    #[allow(dead_code)]
    fn stat(&self) {
        let total_cells = self.rows * self.cols;
        let banded_cells = self.num_cells();
        let percent_cells = (banded_cells as f64) / (total_cells as f64) * 100.0;
        println!(" {} of {} cells are in the band ({2:.2}%)",
                 banded_cells,
                 total_cells,
                 percent_cells);
    }
}


#[cfg(test)]
mod banded {
    use alignment::pairwise::{self, banded, Scoring};
    use alignment::sparse::hash_kmers;
    use utils::TextSlice;

    // Check that the banded alignment is equivalent to the exhaustive SW alignment
    fn compare_to_full_alignment_local(x: TextSlice, y: TextSlice) {

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };


        let mut banded_aligner =
            banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let banded_alignment = banded_aligner.local(x, y);
        // banded_aligner.visualize(&banded_alignment);

        let mut full_aligner = pairwise::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let full_alignment = full_aligner.local(x, y);

        assert_eq!(banded_alignment, full_alignment);
    }

    fn compare_to_full_alignment_global(x: TextSlice, y: TextSlice) {

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };


        let mut banded_aligner =
            banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let banded_alignment = banded_aligner.global(x, y);
        banded_aligner.visualize(&banded_alignment);

        let mut full_aligner = pairwise::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let full_alignment = full_aligner.global(x, y);

        assert_eq!(banded_alignment, full_alignment);
    }

    fn compare_to_full_alignment_semiglobal(x: TextSlice, y: TextSlice) {

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };

        let mut banded_aligner =
            banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let banded_alignment = banded_aligner.semiglobal(x, y);
        banded_aligner.visualize(&banded_alignment);

        let mut full_aligner = pairwise::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score);
        let full_alignment = full_aligner.semiglobal(x, y);
        // banded_aligner.visualize(&full_alignment);

        assert_eq!(banded_alignment, full_alignment);

        let set = hash_kmers(y, 10);
        let banded_alignment = banded_aligner.semiglobal_with_prehash(x, y, &set);
        assert_eq!(banded_alignment, full_alignment);
    }

    #[test]
    fn test_band_add_entry() {
        let mut band = banded::Band::new(10, 10);
        band.add_entry((3u32, 3u32), 3usize);
        assert_eq!(band.ranges,
                   [0..7, 0..7, 0..7, 0..7, 0..7, 0..7, 0..7, 11..0, 11..0, 11..0, 11..0]);
        band.add_entry((9u32, 9u32), 2usize);
        assert_eq!(band.ranges,
                   [0..7, 0..7, 0..7, 0..7, 0..7, 0..7, 0..7, 7..11, 7..11, 7..11, 7..11]);
        band.add_entry((7u32, 5u32), 2usize);
        assert_eq!(band.ranges,
                   [0..7, 0..7, 0..7, 0..10, 0..10, 0..10, 0..10, 5..11, 7..11, 7..11, 7..11]);

        let mut band = banded::Band::new(10, 10);
        band.add_entry((0u32, 0u32), 2usize);
        assert_eq!(band.ranges,
                   [0..3, 0..3, 0..3, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0]);

        let mut band = banded::Band::new(10, 10);
        band.add_entry((10u32, 10u32), 2usize);
        assert_eq!(band.ranges,
                   [11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 8..11, 8..11, 8..11]);

        let mut band = banded::Band::new(10, 10);
        band.add_entry((10u32, 0u32), 2usize);
        assert_eq!(band.ranges,
                   [8..11, 8..11, 8..11, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0]);

        let mut band = banded::Band::new(10, 10);
        band.add_entry((0u32, 10u32), 2usize);
        assert_eq!(band.ranges,
                   [11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 11..0, 0..3, 0..3, 0..3]);
        band.stat();

    }

    fn compare_add_kmer_and_add_entry(start: (u32, u32), k: usize, w: usize, m: usize, n: usize) {
        let mut band1 = banded::Band::new(m, n);
        band1.add_kmer(start, k, w);

        let mut band2 = banded::Band::new(m, n);
        for i in 0..k {
            band2.add_entry((start.0 + i as u32, start.1 + i as u32), w);
        }
        assert_eq!(band1.ranges, band2.ranges);
    }

    #[test]
    fn test_band_add_kmer() {
        compare_add_kmer_and_add_entry((3u32, 3u32), 4, 2, 10, 10);
        compare_add_kmer_and_add_entry((3u32, 3u32), 8, 2, 10, 10);
        compare_add_kmer_and_add_entry((5u32, 0u32), 6, 3, 10, 10);
    }

    #[test]
    fn test_same() {
        let x = b"ACGTATCATAGACCCTAGATAGGGTTGTGTAGATGATCCACAGACGTATCATAGATTAGATAGGGTTGTGTAGATGATTCC\
        ACAG";
        let y = x.clone();
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_big() {
        let query = b"CATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCGTGAGTGACAGCCCCAAGGACACCAAGGGATGAAGCTT\
        CTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTCAGCAGG\
        GCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGGCTCCAT\
        CCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTG";

        let target =
            b"CCTCCCATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCATGAGTGACAGCCCCAAGGACACCAAGGGATG\
        AAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTC\
        AGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGG\
        CTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCG\
        TCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGCACGGCTCCCAA\
        CTCTCTTCCGGCCAAGGATCCCGTGTTCCTGAAATGTCTTTCTACCAAACACAGTTGCTGTGTAACCACTCATTTCATTTTCCTAATTTGT\
        GTTGATCCAGGACACGGGAGGAGACCTGGGCAGCGGCGGACTCATTGCAGGTCGCTCTGCGGTGAGGACGCCACAGGCAC";

        compare_to_full_alignment_local(query, target);
        // compare_to_full_alignment_global(query, target);
        // compare_to_full_alignment_semiglobal(query, target);
    }

    #[test]
    fn test_deletion() {
        let x = b"AGCACACGTGTGCGCTATACAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACACGTGTGCGCTATACAGTAAAAAAAACACGTGTCACAGTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_insertion() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACACGTGTGCGCTATACAGTACACGTGTCACAGTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_substitutions() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_overhangs1() {
        let x = b"CGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"AGCACAAGTGTGAGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_overhangs2() {
        let x = b"GCACACGAGCACACGTAGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGAC";
        let y = b"TATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_overhangs3() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTG";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCACATTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_overhangs4() {
        let x = b"AGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCACAGTTGTACTAGCATGACCAGTTGTACTAGCATGAC";
        let y = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCA";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_overhangs5() {
        let x = b"AGCACAAGTGTGCGCTATACAGGAAGTAGGAGTACACGTGTCA";
        let y = b"CAGTTGTACTAGCATGACCAGTTGTACTAGCATGACAGCACACGTGTGCGCTATACAGTAAGTAGTAGTACACGTGTCA\
            CAGTTGTACTAGCATGACCAGTTGTACTAGCATGAC";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_band_starts_inside() {
        let x = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG";
        let y = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGG";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_band_ends_inside() {
        let x = b"GGGGGGGGGGGGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let y = b"GGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    #[test]
    fn test_band_is_fully_inside() {
        let x = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGG\
        AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let y = b"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTGGGGGGGGGGGGGGGGGGGG\
        TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
        compare_to_full_alignment_local(x, y);
        compare_to_full_alignment_global(x, y);
        compare_to_full_alignment_semiglobal(x, y);
    }

    // #[test]
    // fn test_failure() {
    //     let x = b"AGAATTTTAGTGATCATATCGTTAACAGCTCCTGAGGGGACTTGGCCCAGCTGGAATAGCTACATGGCAGATTTTTCGTTCACGTTTTCTCTTCGCGATGTTCATCACTGCTTCTACCATATCGAGACCCGTTTATTGACTTCAGACAATGAGGAGCAATTAGGACGTTTATACGATGTTAGCGCGTTTAATAACTCACTGATATGCCACAGGCGCAGGCCTGACAAAGTTTATCCGGGGTCGGGAAAGCTGTGCCCTCATCCAAGTGCTCAGCTAACCAGCAACTGTCGGCTAATTCTTAGATATACCGGATTTATTAACACTGGCCTGACATCCTATACCGAGTAGGCCCCCAAAGTAATTGATGTTCCCGCAACTACTACTCCCGAGGCTAGGTCGAGTCCTACTCCAAGACATCCTGCGTAAAGACAAGGCGCTGACTTGACGTAGTAAAGACCTGGCGCGGGATACACACAGCATAGCGTGAAGCACAGACAAACTGAAGTGGCCGAAGAGAATCTAACAATGGTAC";
    //     let y = b"GTTTCGATGCTCACTGAACAGTAGAGTTTACGCCCAACGGTTAGTACCTCGCTAAGGGAGTGGGTGTCCGGGCAGAATTTTAGTGATCATATCGTTAACAGCTCCTGAGGGGACTTGGCCCAGCTGGAATAGCTACATGGCAGATTTTTCGTTCACGTTTTCTCTTCCCGATGTTCATCACTGCTTCTACCATATCGCATCCAAGTGCTCAGCTAACCAGCAACTGTCGGCTAATTCTTAGATATACCGGATTTATTAACACTGGCCTGACATCCTATACCGAGTAGGCCCCCAAAGTAATTGATGTTCCCGCAACTACTACTCCCGAGGCTAGGTATTTGTACCTGTTGCCGCCACGTATCGGGGGCGCTACGGGCGGCACGGCCCGATGCCTTGCTTCCCAGGGTGTTTTTTAGGATTCGATTCAGTGGTCGGTCGGGCTTTAAGCGGTCCAGATCTTAGCTGTATCTCGAGTCCTACTCCAAGACGTCCTGCGTAAAGACAAGGCGCTGACTTGACGTAGTAAAGACCTGGCGCGGGATACACACAGCATAGCGTGAAGCACAGACAAACTGAAGTGGCCGAAGAGAATCTAACAATGGTACTGACAGG";
    //     compare_to_full_alignment_semiglobal(x, y);
    // }

    use alignment::AlignmentOperation::*;
    use std::iter::repeat;
    use scores::blosum62;

    #[test]
    fn test_semiglobal() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }


    // Test case for underflow of the SW score.
    #[test]
    fn test_semiglobal_gap_open_lt_mismatch() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -5i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -1, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Del, Match, Ins, Match, Match, Match]);
    }


    #[test]
    fn test_global_affine_ins() {
        let x = b"ACGAGAACA";
        let y = b"ACGACA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);

        println!("aln:\n{}", alignment.pretty(x, y));
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Ins, Ins, Ins, Match, Match, Match]);
    }

    #[test]
    fn test_global_affine_ins2() {
        let x = b"AGATAGATAGATAGGGAGTTGTGTAGATGATCCACAGT";
        let y = b"AGATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);

        println!("aln:\n{}", alignment.pretty(x, y));

        let mut correct = Vec::new();
        correct.extend(repeat(Match).take(11));
        correct.extend(repeat(Ins).take(10));
        correct.extend(repeat(Match).take(17));

        assert_eq!(alignment.operations, correct);
    }


    #[test]
    fn test_local_affine_ins2() {
        let x = b"ACGTATCATAGATAGATAGGGTTGTGTAGATGATCCACAG";
        let y = b"CGTATCATAGATAGATGTAGATGATCCACAGT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.xstart, 1);
        assert_eq!(alignment.ystart, 0);
    }


    #[test]
    fn test_local() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);
    }

    #[test]
    fn test_global() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);

        println!("\naln:\n{}", alignment.pretty(x, y));
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match,
                    Match]);
    }

    #[test]
    fn test_blosum62() {
        let x = b"AAAA";
        let y = b"AAAA";
        let score = &blosum62;
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, score, 10, 10);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.score, 16);
        assert_eq!(alignment.operations, [Match, Match, Match, Match]);
    }

    #[test]
    fn test_issue11() {
        let y = b"TACC"; //GTGGAC";
        let x = b"AAAAACC"; //GTTGACGCAA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Ins, Ins, Ins, Subst, Match, Match, Match]);
    }


    #[test]
    fn test_issue12_1() {
        let x = b"CCGGCA";
        let y = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 1);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Subst, Subst, Subst]);
    }

    #[test]
    fn test_issue12_2() {
        let y = b"CCGGCA";
        let x = b"ACCGTTGACGC";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.ystart, 0);

        assert_eq!(alignment.operations,
                   [Subst, Match, Ins, Ins, Ins, Ins, Ins, Ins, Subst, Match, Match]);
    }


    #[test]
    fn test_issue12_3() {
        let y = b"CCGTCCGGCAA";
        let x = b"AAAAACCGTTGACGCAA";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Ins, Ins, Ins, Ins, Ins, Ins, Match, Subst, Subst, Match, Subst, Subst,
                    Subst, Match, Match, Match, Match]);


        let mut aligner = banded::Aligner::with_capacity(y.len(), x.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(y, x);

        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Subst, Subst, Match, Subst, Subst, Subst, Match, Match, Match, Match]);
    }


    #[test]
    fn test_left_aligned_del() {
        let x = b"GTGCATCATGTG";
        let y = b"GTGCATCATCATGTG";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Del, Del, Del, Match, Match, Match, Match, Match, Match,
                    Match, Match, Match]);
    }


    // Test that trailing deletions are correctly handled
    // in global mode
    #[test]
    fn test_global_right_del() {
        let x = b"AACCACGTACGTGGGGGGA";
        let y = b"CCACGTACGT";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);

        println!("\naln:\n{}", alignment.pretty(x, y));

        println!("score:{}", alignment.score);
        assert_eq!(alignment.score, -9);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Ins, Ins, Match, Match, Match, Match, Match, Match, Match, Match, Match,
                    Match, Ins, Ins, Ins, Ins, Ins, Ins, Ins]);
    }



    #[test]
    fn test_left_aligned_ins() {

        let x = b"GTGCATCATCATGTG";
        let y = b"GTGCATCATGTG";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::with_capacity(x.len(), y.len(), -5, -1, &score, 10, 10);
        let alignment = aligner.global(x, y);
        println!("\naln:\n{}", alignment.pretty(x, y));

        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Ins, Ins, Ins, Match, Match, Match, Match, Match, Match,
                    Match, Match, Match]);
    }



    #[test]
    fn test_aligner_new() {
        let x = b"ACCGTGGAT";
        let y = b"AAAAACCGTTGAT";
        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::new(-5, -1, &score, 10, 10);

        let alignment = aligner.semiglobal(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);

        let alignment = aligner.local(x, y);
        assert_eq!(alignment.ystart, 4);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Match, Match, Match, Match, Match, Subst, Match, Match, Match]);

        let alignment = aligner.global(x, y);
        assert_eq!(alignment.ystart, 0);
        assert_eq!(alignment.xstart, 0);
        assert_eq!(alignment.operations,
                   [Del, Del, Del, Del, Match, Match, Match, Match, Match, Subst, Match, Match,
                    Match]);
    }

    #[test]
    fn test_semiglobal_simple() {

        let x = b"GAAAACCGTTGAT";
        let y = b"ACCGTGGATGGG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let mut aligner = banded::Aligner::new(-5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(alignment.operations,
                   [Ins, Ins, Ins, Ins, Match, Match, Match, Match, Match, Subst, Match, Match,
                    Match]);
    }

    #[test]
    fn test_insert_only_semiglobal() {

        let x = b"TTTT";
        let y = b"AAAA";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let mut aligner = banded::Aligner::new(-5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(alignment.operations, [Ins, Ins, Ins, Ins]);
    }

    #[test]
    fn test_insert_in_between_semiglobal() {

        let x = b"GGGGG";
        let y = b"GGTAGGG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let mut aligner = banded::Aligner::new(-5, -1, &score, 10, 10);
        let alignment = aligner.semiglobal(x, y);

        assert_eq!(alignment.operations,
                   [Match, Match, Del, Del, Match, Match, Match]);
    }

    #[test]
    fn test_xclip_prefix_custom() {

        let x = b"GGGGGGATG";
        let y = b"ATG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5);

        let mut aligner = banded::Aligner::with_scoring(scoring, 10, 10);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Xclip(6), Match, Match, Match]);
    }

    #[test]
    fn test_yclip_prefix_custom() {

        let y = b"GGGGGGATG";
        let x = b"ATG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).yclip(-5);

        let mut aligner = banded::Aligner::with_scoring(scoring, 10, 10);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Yclip(6), Match, Match, Match]);
    }

    #[test]
    fn test_xclip_suffix_custom() {

        let x = b"GAAAA";
        let y = b"CG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -1i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5).yclip(0);

        let mut aligner = banded::Aligner::with_scoring(scoring, 10, 10);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Yclip(1), Match, Xclip(4)]);
    }

    #[test]
    fn test_yclip_suffix_custom() {

        let y = b"GAAAA";
        let x = b"CG";

        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let scoring = Scoring::new(-5, -1, &score).yclip(-5).xclip(0);

        let mut aligner = banded::Aligner::with_scoring(scoring, 10, 10);
        let alignment = aligner.custom(x, y);

        assert_eq!(alignment.operations, [Xclip(1), Match, Yclip(4)]);
    }

    #[test]
    fn test_longer_string_all_operations() {

        let x = b"TTTTTGGGGGGATGGCCCCCCTTTTTTTTTTGGGAAAAAAAAAGGGGGG";
        let y = b"GGGGGGATTTCCCCCCCCCTTTTTTTTTTAAAAAAAAA";


        let score = |a: u8, b: u8| if a == b { 1i32 } else { -3i32 };
        let scoring = Scoring::new(-5, -1, &score).xclip(-5).yclip(0);

        let mut aligner = banded::Aligner::with_scoring(scoring, 10, 10);
        let alignment = aligner.custom(x, y);

        println!("{}", alignment.pretty(x, y));
        assert_eq!(alignment.score, 7);

    }
}

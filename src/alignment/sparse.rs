// Copyright 2014-2015 Patrick Marks
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Calculate 'sparse' alignments from kmer matches. Can be much faster than
//! Smith-Waterman for long string, when a large enough k is used.
//! Complexity: O(n * log(n)) for a pair of strings with n k-kmer matches. This
//! approach is useful for generating an approximate 'backbone' alignments
//! between two long sequences, for example in long-read alignment or
//! genome-genome alignment. The backbone alignment can be used as-is, or can serve
//! as a guide for a banded alignment.  By tuning k so that len(query) + len(reference) < 4^k,
//! the number of false positive kmer matches is kept small, resulting in very
//! fast run times for long strings.
//!
//! # Example
//!
//! ```
//! use bio::alignment::sparse::*;
//!
//! let s1 =   b"ACGTACGATAGGTA";
//! let s2 = b"TTACGTACGATAGGTATT";
//! let k = 8;
//! let matches = find_kmer_matches(s1, s2, k);
//! let sparse_al = lcskpp(&matches, k);
//! let match_path: Vec<(u32,u32)> = sparse_al.path.iter().map(|i| matches[*i]).collect();
//! assert_eq!(match_path, vec![(0,2), (1,3), (2,4), (3,5), (4,6), (5,7), (6,8)]);
//! assert_eq!(sparse_al.score, 14);

use crate::data_structures::bit_tree::MaxBitTree;
use fxhash::FxHasher;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::hash::BuildHasherDefault;

pub type HashMapFx<K, V> = HashMap<K, V, BuildHasherDefault<FxHasher>>;

/// Result of a sparse alignment
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct SparseAlignmentResult {
    /// LCSk++ path, represented as vector of indices into the input matches vector.
    pub path: Vec<usize>,
    // Score of the path, which is the number of bases covered by the matched kmers.
    pub score: u32,
    // Full DP vector, which can generally be ignored. (It may be useful for testing purposes).
    pub dp_vector: Vec<(u32, i32)>,
}

/// Sparse DP routine for Longest Common Subsequence in length k substrings.  Also known of LCSk++
/// From LCSk++: Practical similarity metric for long strings. Filip Pavetić, Goran Žužić, Mile Šikić
/// Paper here :https://arxiv.org/abs/1407.2407.  Original implementation here:
/// https://github.com/fpavetic/lcskpp
///
/// # Arguments
///
/// * `matches` - a vector of tuples indicating the (string1 position, string2 position) kmer
///   matches between the strings
/// * `k` - the kmer length used for matching
///
/// # Return value
///
/// The method returns a `SparseAlignmentResult` struct with the following fields:
/// * `path` is the LCSk++ path, represented as vector of indices into the input matches vector.
/// * `score` is the score of the path, which is the number of bases covered by the matched kmers.
/// * `dp_vector` is the full DP vector, which can generally be ignored. (It may be useful for
///   testing purposes).
pub fn lcskpp(matches: &[(u32, u32)], k: usize) -> SparseAlignmentResult {
    if matches.is_empty() {
        return SparseAlignmentResult {
            path: Vec::new(),
            score: 0,
            dp_vector: Vec::new(),
        };
    }

    let k = k as u32;

    // incoming matches must be sorted to let us find the predecessor kmers by binary search.
    for i in 1..matches.len() {
        assert!(matches[i - 1] < matches[i]);
    }

    let mut events: Vec<(u32, u32, u32)> = Vec::new();
    let mut n = 0;

    for (idx, &(x, y)) in matches.iter().enumerate() {
        events.push((x, y, (idx + matches.len()) as u32));
        events.push((x + k, y + k, idx as u32));

        n = max(n, x + k);
        n = max(n, y + k);
    }
    events.sort_unstable();

    let mut max_col_dp: MaxBitTree<(u32, u32)> = MaxBitTree::new(n as usize);
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0);

    dp.resize(events.len(), (0, 0));

    for ev in events {
        let p = (ev.2 % matches.len() as u32) as usize;
        let j = ev.1;
        let is_start = ev.2 >= (matches.len() as u32);

        if is_start {
            dp[p] = (k, -1);
            let (best_value, best_position) = max_col_dp.get(j as usize);
            if best_value > 0 {
                dp[p] = (k + best_value, best_position as i32);
                best_dp = max(best_dp, (dp[p].0, p as i32));
            }
        } else {
            // See if this kmer continues a different kmer
            if ev.0 > k && ev.1 > k {
                if let Ok(cont_idx) = matches.binary_search(&(ev.0 - k - 1, ev.1 - k - 1)) {
                    let prev_score = dp[cont_idx].0;
                    let candidate = (prev_score + 1, cont_idx as i32);
                    dp[p] = max(dp[p], candidate);
                    best_dp = max(best_dp, (dp[p].0, p as i32));
                }
            }

            max_col_dp.set(ev.1 as usize, (dp[p].0, p as u32));
        }
    }

    let mut traceback = Vec::new();
    let (best_score, mut prev_match) = best_dp;
    while prev_match >= 0 {
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    traceback.reverse();
    SparseAlignmentResult {
        path: traceback,
        score: best_score,
        dp_vector: dp,
    }
}

#[derive(PartialEq, Eq, Ord, PartialOrd, Default, Copy, Clone)]
struct PrevPtr {
    plane: u32,
    score: u32,
    d: u32,
    id: usize,
    x: u32,
    y: u32,
}

impl PrevPtr {
    pub fn new(score: u32, x: u32, y: u32, id: usize, gap_extend: u32) -> PrevPtr {
        let d = x + y;
        PrevPtr {
            plane: score + (d * gap_extend),
            score,
            d,
            id,
            x,
            y,
        }
    }
}

/// Sparse DP routine generalizing LCSk++ method above to penalize alignment gaps.
/// A gap is an unknown combination of mismatch, insertion and deletions, and incurs
/// a penalty of gap_open + d * gap_extend, where d is the distance along the diagonal of the gap.
/// # Arguments
///
/// * `matches` - a vector of tuples indicating the (string1 position, string2 position) kmer
///   matches between the strings
/// * `k` - the kmer length used for matching
/// * `match_score` - reward for each matched base
/// * `gap_open` - score of opening a gap, including a mismatch gap. Must be negative.
/// * `gap_extend` - score for extending a gap along the diagonal. Must be negative.
///
/// # Return value
///
/// The method returns a `SparseAlignmentResult` struct with the following fields:
/// * `path` is the SDP path, represented as vector of indices into the input matches vector.
/// * `score` is the score of the path, which is the number of bases covered by the matched kmers.
/// * `dp_vector` is the full DP vector, which can generally be ignored. (It may be useful for
///   testing purposes).
pub fn sdpkpp(
    matches: &[(u32, u32)],
    k: usize,
    match_score: u32,
    gap_open: i32,
    gap_extend: i32,
) -> SparseAlignmentResult {
    if matches.is_empty() {
        return SparseAlignmentResult {
            path: Vec::new(),
            score: 0,
            dp_vector: Vec::new(),
        };
    }

    let k = k as u32;
    if gap_open > 0 || gap_extend > 0 {
        panic!("gap parameters cannot be positive")
    }
    let _gap_open = (-gap_open) as u32;
    let _gap_extend = (-gap_extend) as u32;

    // incoming matches must be sorted to let us find the predecessor kmers by binary search.
    for i in 1..matches.len() {
        assert!(matches[i - 1] < matches[i]);
    }

    let mut events: Vec<(u32, u32, u32)> = Vec::new();
    let mut n = 0;

    for (idx, &(x, y)) in matches.iter().enumerate() {
        events.push((x, y, (idx + matches.len()) as u32));
        events.push((x + k, y + k, idx as u32));

        n = max(n, x + k);
        n = max(n, y + k);
    }
    events.sort_unstable();

    let mut max_col_dp: MaxBitTree<PrevPtr> = MaxBitTree::new(n as usize);
    let mut dp: Vec<(u32, i32)> = Vec::with_capacity(events.len());
    let mut best_dp = (k, 0);

    dp.resize(events.len(), (0, 0));

    for ev in events {
        let p = (ev.2 % matches.len() as u32) as usize;
        let j = ev.1;
        let is_start = ev.2 >= (matches.len() as u32);

        if is_start {
            // Default case -- chain starts at this node
            dp[p] = (k * match_score, -1);

            // Find best previous chain, and extend.
            let best_prev = max_col_dp.get(j as usize);
            if best_prev.score > 0 {
                let prev_x = best_prev.x;
                let prev_y = best_prev.y;
                let cur_x = ev.0;
                let cur_y = ev.1;
                let gap = max(cur_x - prev_x, cur_y - prev_y);
                let gap_penalty = if gap > 0 {
                    _gap_open + gap * _gap_extend
                } else {
                    0
                };

                let reward = k * match_score;
                let new_score = (best_prev.score + reward).saturating_sub(gap_penalty);

                dp[p] = max(dp[p], (new_score, best_prev.id as i32));
                best_dp = max(best_dp, (dp[p].0, p as i32));
            }
        } else {
            // See if this kmer continues a different kmer
            if ev.0 > k && ev.1 > k {
                if let Ok(cont_idx) = matches.binary_search(&(ev.0 - k - 1, ev.1 - k - 1)) {
                    let prev_score = dp[cont_idx].0;
                    let candidate = (prev_score + match_score, cont_idx as i32);
                    dp[p] = max(dp[p], candidate);
                    best_dp = max(best_dp, (dp[p].0, p as i32));
                }
            }

            let prev_frag = PrevPtr::new(dp[p].0, ev.0, ev.1, p, _gap_extend);
            max_col_dp.set(ev.1 as usize, prev_frag);
        }
    }

    let mut traceback = Vec::new();
    let (best_score, mut prev_match) = best_dp;
    while prev_match >= 0 {
        traceback.push(prev_match as usize);
        prev_match = dp[prev_match as usize].1;
    }
    traceback.reverse();
    SparseAlignmentResult {
        path: traceback,
        score: best_score,
        dp_vector: dp,
    }
}

pub fn sdpkpp_union_lcskpp_path(
    matches: &[(u32, u32)],
    k: usize,
    match_score: u32,
    gap_open: i32,
    gap_extend: i32,
) -> Vec<usize> {
    if matches.is_empty() {
        return Vec::new();
    }
    let lcskpp_al = lcskpp(matches, k);
    let sdpkpp_al = sdpkpp(matches, k, match_score, gap_open, gap_extend);
    let pre_lcskpp = lcskpp_al
        .path
        .binary_search(&sdpkpp_al.path[0])
        .unwrap_or(0);
    let post_lcskpp = match lcskpp_al.path.binary_search(sdpkpp_al.path.last().unwrap()) {
        Ok(ind) => ind + 1,
        Err(_) => lcskpp_al.path.len(),
    };

    let mut path_union = Vec::new();
    for i in 0..pre_lcskpp {
        path_union.push(lcskpp_al.path[i]);
    }
    for i in 0..sdpkpp_al.path.len() {
        path_union.push(sdpkpp_al.path[i]);
    }
    for i in post_lcskpp..lcskpp_al.path.len() {
        path_union.push(lcskpp_al.path[i]);
    }

    path_union
}

/// Find all matches of length k between two strings, using a q-gram
/// index. For very long reference strings, it may be more efficient to use and
/// FMD index to generate the matches. Note that this method is mainly for
/// demonstration & testing purposes.  For aligning many query sequences
/// against the same reference, you should reuse the QGramIndex of the reference.
pub fn find_kmer_matches(seq1: &[u8], seq2: &[u8], k: usize) -> Vec<(u32, u32)> {
    if seq1.len() < seq2.len() {
        let set = hash_kmers(seq1, k);
        find_kmer_matches_seq1_hashed(&set, seq2, k)
    } else {
        let set = hash_kmers(seq2, k);
        find_kmer_matches_seq2_hashed(seq1, &set, k)
    }
}

/// Creates a HashMap containing all the k-mers in the sequence. FxHasher is used
/// as the hash function instead of the inbuilt one. A good rolling hash function
/// should speed up the code.
pub fn hash_kmers(seq: &[u8], k: usize) -> HashMapFx<&[u8], Vec<u32>> {
    let slc = seq;
    let mut set: HashMapFx<&[u8], Vec<u32>> = HashMapFx::default();
    for i in 0..(slc.len() + 1).saturating_sub(k) {
        set.entry(&slc[i..i + k]).or_default().push(i as u32);
    }
    set
}

// Find all matches of length k between two strings where the first string is
// already hashed by using the function sparse::hash_kmers
pub fn find_kmer_matches_seq1_hashed(
    seq1_set: &HashMapFx<&[u8], Vec<u32>>,
    seq2: &[u8],
    k: usize,
) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();

    for i in 0..(seq2.len() + 1).saturating_sub(k) {
        let slc = &seq2[i..i + k];
        if let Some(matches1) = seq1_set.get(slc) {
            for pos1 in matches1 {
                matches.push((*pos1, i as u32));
            }
        }
    }

    matches.sort_unstable();
    matches
}

// Find all matches of length k between two strings where the second string is
// already hashed by using the function sparse::hash_kmers
pub fn find_kmer_matches_seq2_hashed(
    seq1: &[u8],
    seq2_set: &HashMapFx<&[u8], Vec<u32>>,
    k: usize,
) -> Vec<(u32, u32)> {
    let mut matches = Vec::new();

    for i in 0..(seq1.len() + 1).saturating_sub(k) {
        let slc = &seq1[i..i + k];

        if let Some(matches1) = seq2_set.get(slc) {
            for pos1 in matches1 {
                matches.push((i as u32, *pos1));
            }
        }
    }

    matches.sort_unstable();
    matches
}

// Search for mismatches between two sequences, using an array of sorted matches as input.
// First proceed with expanding by 1 unit on the left side of the matches,
// then proceed on an right expand of k+1 unit
// A maximum number of mismatches can be input, the test is a strict one, which mean for 
// allowed_mismatches = 1, a maximum of 2 mismatches can be found before stopping.
pub fn expand_kmer_matches(
    seq1: &[u8],
    seq2: &[u8],
    k: usize,
    sorted_matches: &[(u32, u32)],
    allowed_mismatches: usize,
) -> Vec<(u32, u32)> {
    // incoming matches must be sorted.
    for i in 1..sorted_matches.len() {
        assert!(sorted_matches[i - 1] < sorted_matches[i]);
    }

    let mut last_match_along_diagonal: HashMapFx<i32, (i32, i32)> = HashMapFx::default();
    let mut left_expanded_matches: Vec<(u32, u32)> = sorted_matches.to_owned();

    for &this_match in sorted_matches.iter() {
        let diag = (this_match.0 as i32) - (this_match.1 as i32);
        let min_xy = min(this_match.0, this_match.1) as i32;
        let default_last_match = (
            this_match.0 as i32 - min_xy - 1,
            this_match.1 as i32 - min_xy - 1,
        );
        let last_match = last_match_along_diagonal
            .get(&diag)
            .cloned()
            .unwrap_or(default_last_match);

        let mut n_mismatches = 0;
        let mut curr_pos = (this_match.0 as i32 - 1, this_match.1 as i32 - 1);
        loop {
            // Reached last position of current match
            if last_match >= curr_pos {
                break;
            }
            // Add the current position mismatch to the total sum of mismatch
            n_mismatches += if seq1[curr_pos.0 as usize] == seq2[curr_pos.1 as usize] {
                0
            } else {
                1
            };
            // Maximum number of mismatches as been reached
            if n_mismatches > allowed_mismatches {
                break;
            }
            left_expanded_matches.push((curr_pos.0 as u32, curr_pos.1 as u32));
            curr_pos = (curr_pos.0 - 1, curr_pos.1 - 1);
        }
        // We need to check until 1 position after this match, when we start our search from
        // the next kmer match along this diagonal
        last_match_along_diagonal.insert(diag, (this_match.0 as i32, this_match.1 as i32));
    }

    left_expanded_matches.sort_unstable();
    let mut expanded_matches = left_expanded_matches.clone();
    left_expanded_matches.reverse();

    let mut next_match_along_diagonal: HashMapFx<i32, (u32, u32)> = HashMapFx::default();

    // Reverting the left expanded_matches to test on right expand
    for &this_match in &left_expanded_matches {
        let diag = (this_match.0 as i32) - (this_match.1 as i32);
        let max_inc = (min(
            seq1.len() as u32 - this_match.0,
            seq2.len() as u32 - this_match.1,
        ) as u32)
            .saturating_sub(k as u32 - 1);
        let next_match = next_match_along_diagonal
            .get(&diag)
            .cloned()
            .unwrap_or((this_match.0 + max_inc, this_match.1 + max_inc));

        let mut n_mismatches = 0;
        let mut curr_pos = (this_match.0 + 1, this_match.1 + 1);
        loop {
            // println!(" This : ({},{}), Current : ({},{}), Next : ({}, {}), Miss : {}",
            // this_match.0, this_match.1, curr_pos.0, curr_pos.1, next_match.0, next_match.1, n_mismatches);
            if curr_pos >= next_match {
                break;
            }
            n_mismatches +=
                if seq1[curr_pos.0 as usize + k - 1] == seq2[curr_pos.1 as usize + k - 1] {
                    0
                } else {
                    1
                };
            if n_mismatches > allowed_mismatches {
                break;
            }
            expanded_matches.push(curr_pos);
            curr_pos = (curr_pos.0 + 1, curr_pos.1 + 1);
        }

        next_match_along_diagonal.insert(diag, this_match);
    }
    expanded_matches.sort_unstable();
    expanded_matches
}

#[cfg(test)]
mod sparse_alignment {
    use super::find_kmer_matches;

    #[test]
    pub fn test_find_kmer_matches() {
        let s1 = b"ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = b"TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        let k = 8;
        //let s1 = "  ACGTACGATAGATCCGTACGTAACA     GTACAGTATATCAGTTATATGCGATA";
        //let s2 = "TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";

        let hits = find_kmer_matches(s1, s2, k);
        assert_eq!(hits.len(), (25 - k + 1) + (24 - k + 1));
        //println!("hits: {:?}", hits);
    }

    #[test]
    pub fn test_lcskpp0() {
        let s1 = b"ACGTACGATAGGTA";
        let s2 = b"TTACGTACGATAGGTATT";
        let k = 8;
        let matches = super::find_kmer_matches(s1, s2, k);
        let res = super::lcskpp(&matches, k);
        let match_path: Vec<(u32, u32)> = res.path.iter().map(|i| matches[*i]).collect();
        assert_eq!(
            match_path,
            vec![(0, 2), (1, 3), (2, 4), (3, 5), (4, 6), (5, 7), (6, 8)]
        );
        assert_eq!(res.score, 14);
    }

    pub fn strict_compare_lcskpp_sdpkpp(s1: &[u8], s2: &[u8]) {
        let k = 8;
        let matches = super::find_kmer_matches(s1, s2, k);
        let res1 = super::lcskpp(&matches, k);
        let res2 = super::sdpkpp(&matches, k, 1, 0, 0);

        assert_eq!(res1, res2);
    }

    #[test]
    pub fn test_sdp() {
        let s1 = b"ACGTACGATAGGTA";
        let s2 = b"TTACGTACGATAGGTATT";
        strict_compare_lcskpp_sdpkpp(s1, s2);
    }

    #[test]
    pub fn test_lcskpp1() {
        let s1 = b"ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = b"TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        let k = 8;

        let matches = super::find_kmer_matches(s1, s2, k);
        let res = super::lcskpp(&matches, k);

        // For debugging:
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        // Should have 25bp group of matches plus a 24 bp group of matches
        assert_eq!(res.score, 25 + 24);
    }

    #[test]
    pub fn test_sdp1() {
        let s1 = b"ACGTACGATAGATCCGTACGTAACAGTACAGTATATCAGTTATATGCGATA";
        let s2 = b"TTACGTACGATAGATCCGTACGTAACATTTTTGTACAGTATATCAGTTATATGCGA";
        strict_compare_lcskpp_sdpkpp(s1, s2);
    }

    #[test]
    pub fn test_lcskpp2() {
        // Match the same string -- should get a diagonal traceback, despite lots of off-diagonal
        // homology
        let s1 = b"ACGTACGATAGATCCGACGTACGTACGTTCAGTTATATGACGTACGTACGTAACATTTTTGTA";
        let k = 5;

        let matches = super::find_kmer_matches(s1, s1, k);
        let res = super::lcskpp(&matches, k);

        // For debugging:
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        assert_eq!(res.score, s1.len() as u32);

        for i in 0..res.path.len() {
            assert_eq!(matches[res.path[i] as usize], (i as u32, i as u32));
        }
    }

    #[test]
    pub fn test_sdp2() {
        let s1 = b"ACGTACGATAGATCCGACGTACGTACGTTCAGTTATATGACGTACGTACGTAACATTTTTGTA";
        strict_compare_lcskpp_sdpkpp(s1, s1);
    }

    // Test case from local SV caller alignments.
    // The query sequence ends in 1-2 copies of tandem repeat element
    // The target sequence end in >2 copies of the element
    // Without a gap penalty (i.e. LCSk++), the alignment of the
    // TRs is arbitrary, and way the implementation breaks ties may introduce
    // a gap while maintaining the same score.
    // The SDP code with gap open & extend penalties should resolve this.
    const QUERY_REPEAT: &[u8] = b"CCTCCCATCTCCACCCACCCTATCCAACCCTGGGGTGGCAGGTCATGAGTGA\
CAGCCCCAAGGACACCAAGGGATGAAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGCAGAACAGGAGGCCCCATCTCC\
CGTTCTTACTCAGAAGCTGTCAGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTCCCCAGCTCTTGCTTCAGGGCCTCACGTTTCG\
CCCCCTGAGGCCTGGGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTG";

    const TARGET_REPEAT: &[u8] = b"CCTCCCATCTCCACCCACCCTATCCAACCCTGGGGTGGCAG\
GTCATGAGTGACAGCCCCAAGGACACCAAGGGATGAAGCTTCTCCTGTGCTGAGATCCTTCTCGGACTTTCTGAGAGGCCACGC\
AGAACAGGAGGCCCCATCTCCCGTTCTTACTCAGAAGCTGTCAGCAGGGCTGGGCTCAAGATGAACCCGTGGCCGGCCCCACTC\
CCCAGCTCTTGCTTCAGGGCCTCACGTTTCGCCCCCTGAGGCCTGGGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACAT\
CTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAAC\
ATCTGGTGGGCTCCGTCCTCACGGCTGGAGGGGCTCTCAGAACATCTGGTGCACGGCTCCCAACTCTCTTCCGGCCAAGGATCC\
CGTGTTCCTGAAATGTCTTTCTACCAAACACAGTTGCTGTGTAACCACTCATTTCATTTTCCTAATTTGTGTTGATCCAGGACA\
CGGGAGGAGACCTGGGCAGCGGCGGACTCATTGCAGGTCGCTCTGCGGTGAGGACGCCACAGGCAC";

    #[test]
    fn test_lcskpp_tandem_repeat() {
        let k = 8;
        let matches = super::find_kmer_matches(QUERY_REPEAT, TARGET_REPEAT, k);
        let res = super::lcskpp(&matches, k);

        // For debugging:
        //for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
        //    println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        //}
        //println!("tb: {:?}", tb);

        assert_eq!(res.score, QUERY_REPEAT.len() as u32);

        // NOTE -- this test will fail, because LCSk++ introduces a gap in the placement of the TR
        // Corrected with gap scoring in SDP
        /*
        for i in 0..res.path.len() {
            assert_eq!(matches[res.path[i] as usize], (i as u32, i as u32));
        }
        */
    }

    #[test]
    fn test_sdpkpp_tandem_repeat() {
        let k = 8;
        let matches = super::find_kmer_matches(QUERY_REPEAT, TARGET_REPEAT, k);
        let res = super::sdpkpp(&matches, k, 1, -1, -1);

        // For debugging:
        /*
        for (idx, (ev, (score, prev))) in evs.iter().zip(dps.clone()).enumerate() {
            println!("idx: {:?}\tev: {:?}\tscore: {:?}\t prev: {:?}", idx, ev, score, prev);
        }
        println!("tb: {:?}", tb);
        */

        assert_eq!(res.score, QUERY_REPEAT.len() as u32);

        for i in 0..res.path.len() {
            assert_eq!(matches[res.path[i] as usize], (i as u32, i as u32));
        }
    }

    #[test]
    fn test_sdpkpp_same() {
        let x = b"ACGTACGTAC";
        let y = b"ACGTACGTAC";
        let matches = super::find_kmer_matches(x, y, 10);
        let res = super::sdpkpp(&matches, 10, 1, -1, -1);
        assert_eq!(res.path, [0]);
        assert_eq!(res.score, 10);

        let x = b"ACGTACGTACA";
        let y = b"ACGTACGTACA";
        let matches = super::find_kmer_matches(x, y, 10);
        let res = super::sdpkpp(&matches, 10, 1, -1, -1);
        assert_eq!(res.path, [0, 1]);
        assert_eq!(res.score, 11);

        let x = b"ACGTACGTACACGTACGTAC";
        let y = b"ACGTACGTAC";
        let matches = super::find_kmer_matches(x, y, 10);
        let res = super::sdpkpp(&matches, 10, 1, -1, -1);
        assert_eq!(res.path, [0]);
        assert_eq!(res.score, 10);
    }

    #[test]
    fn test_lcskpp_same() {
        let x = b"ACGTACGTAC";
        let y = b"ACGTACGTAC";
        let matches = super::find_kmer_matches(x, y, 10);
        let res = super::lcskpp(&matches, 10);
        assert_eq!(res.path, [0]);
        assert_eq!(res.score, 10);

        let x = b"ACGTACGTACA";
        let y = b"ACGTACGTACA";
        let matches = super::find_kmer_matches(x, y, 10);
        let res = super::lcskpp(&matches, 10);
        assert_eq!(res.path, [0, 1]);
        assert_eq!(res.score, 11);

        let x = b"ACGTACGTACACGTACGTAC";
        let y = b"ACGTACGTAC";
        let matches = super::find_kmer_matches(x, y, 10);
        let res = super::lcskpp(&matches, 10);
        assert_eq!(res.path, [0]);
        assert_eq!(res.score, 10);
    }

    #[test]
    fn test_expanded_matches() {
        let x = b"GGGCAAAAAA";
        let y = b"GGGGAAAAAA";
        let matches = super::find_kmer_matches(x, y, 6);
        assert_eq!(matches, vec![(4, 4)]);

        let expanded_matches = super::expand_kmer_matches(x, y, 6, &matches, 1);
        assert_eq!(
            expanded_matches,
            (0..5).map(|x| (x, x)).collect::<Vec<(u32, u32)>>()
        );

        let x = b"TTTTTTGGGCAAAAAA";
        let y = b"TTTTTTGGGGAAAAAA";
        let matches = super::find_kmer_matches(x, y, 6);
        assert_eq!(matches, vec![(0, 0), (1, 1), (2, 2), (3, 3), (10, 10)]);

        let expanded_matches = super::expand_kmer_matches(x, y, 6, &matches, 1);
        assert_eq!(
            expanded_matches,
            (0..11).map(|x| (x, x)).collect::<Vec<(u32, u32)>>()
        );

        let x = b"TTTTTTCCGCAAAAAA";
        let y = b"TTTTTTGGGGAAAAAA";
        let matches = super::find_kmer_matches(x, y, 6);
        assert_eq!(matches, vec![(0, 0), (10, 10)]);

        let expanded_matches = super::expand_kmer_matches(x, y, 6, &matches, 1);
        assert_eq!(
            expanded_matches,
            vec![(0, 0), (1, 1), (8, 8), (9, 9), (10, 10)]
        );

        let x = b"TTTTTTCGGCAAAAAA";
        let y = b"TTTTTTGGGGAAAAAA";
        let matches = super::find_kmer_matches(x, y, 6);
        assert_eq!(matches, vec![(0, 0), (10, 10)]);

        let expanded_matches = super::expand_kmer_matches(x, y, 6, &matches, 1);
        assert_eq!(
            expanded_matches,
            vec![
                (0, 0),
                (1, 1),
                (2, 2),
                (3, 3),
                (7, 7),
                (8, 8),
                (9, 9),
                (10, 10),
            ]
        );

        let x = b"AAAAAACGGG";
        let y = b"AAAAAAGGGG";
        let matches = super::find_kmer_matches(x, y, 6);
        assert_eq!(matches, vec![(0, 0)]);
        let expanded_matches = super::expand_kmer_matches(x, y, 6, &matches, 1);
        assert_eq!(
            expanded_matches,
            (0..5).map(|x| (x, x)).collect::<Vec<(u32, u32)>>()
        );
    }
}

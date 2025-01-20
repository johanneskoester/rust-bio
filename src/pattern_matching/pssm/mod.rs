// Copyright 2018 Kieran Hervold
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Create a weight matrix representing a set of aligned reference sequences
//! that constitute a motif, and use this matrix to scan query sequences for
//! occurrences of this motif.
//! Complexity: O(n*m) for motif length n and query length m
//!
//! The position-specific scoring matrix (PSSM), aka position weight matrix (PWM),
//! algorithm is implemented for both DNA and amino-acid sequences.
//!
//! # Examples
//!
//! use bio::pattern_matching::pssm::DNAMotif;
//! let pssm = DNAMotif::from_seqs(vec![
//!            b"AAAA".to_vec(),
//!            b"AATA".to_vec(),
//!            b"AAGA".to_vec(),
//!            b"AAAA".to_vec(),
//!        ].as_ref(), None).unwrap();
//! let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
//! println!("motif found at position {}", start_pos);
//!
//! /* amino acid sequences are supported, too */
//! use pssm::pattern_matching::pssm::ProtMotif;
//! let pssm = ProtMotif::from_seqs(vec![
//!            b"ARNNYM".to_vec(),
//!            b"ARNRYM".to_vec(),
//!            b"ARNNCM".to_vec(),
//!            b"ARNNYM".to_vec(),
//!        ].as_ref(), None).unwrap();

use std::borrow::Borrow;
use std::f32::NEG_INFINITY;

use itertools::Itertools;
use ndarray::prelude::*;

mod dnamotif;
pub mod errors;
mod protmotif;

pub use self::dnamotif::DNAMotif;
pub use self::errors::{Error, Result};
pub use self::protmotif::ProtMotif;

/// default pseudocount - used to prevent 0 tallies
pub const DEF_PSEUDO: f32 = 0.5;
/// approximately zero
pub const EPSILON: f32 = 1e-5;
/// value representing an invalid monomer in lookup table
pub const INVALID_MONO: u8 = 255;

/// Represents motif score & location of match
#[derive(Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct ScoredPos {
    pub loc: usize,
    pub sum: f32,
    pub scores: Vec<f32>,
}

impl Default for ScoredPos {
    fn default() -> ScoredPos {
        ScoredPos {
            loc: 0,
            sum: NEG_INFINITY,
            scores: Vec::new(),
        }
    }
}

/// Trait containing code shared between DNA and protein implementations
/// of the position-specific scoring matrix.
pub trait Motif {
    /// Lookup table mapping monomer -> index
    const LK: [u8; 127] = [INVALID_MONO; 127];
    /// All monomers, in order corresponding to lookup table
    const MONOS: &'static [u8] = b"";
    /// Monomer count - equal to length of `MONOS`
    const MONO_CT: usize = 0;

    /// Returns a weight matrix representing the sequences provided.
    /// This code is shared by implementations of `from_seqs`
    /// # Arguments
    /// * `seqs` - sequences incorporated into motif
    /// * `pseudos` - array slice with a pseudocount for each monomer;
    ///    defaults to DEF_PSEUDO for all if None is supplied
    ///
    /// FIXME: pseudos should be an array of size MONO_CT, but that
    /// is currently unsupported in stable rust
    fn seqs_to_weights(seqs: &[Vec<u8>], _pseudos: Option<&[f32]>) -> Result<Array2<f32>> {
        if _pseudos.is_some() {
            if _pseudos.unwrap().len() != Self::MONO_CT {
                return Err(Error::InvalidPseudos {
                    expected: Self::MONO_CT as u8,
                    received: _pseudos.unwrap().len() as u8,
                });
            }
        }
        let pseudos: Array1<f32> = match _pseudos {
            Some(p) => Array::from_vec(p.iter().cloned().collect()),
            None => Array::from_vec(vec![DEF_PSEUDO; Self::MONO_CT]),
        };

        if seqs.is_empty() {
            return Err(Error::EmptyMotif);
        }

        let seqlen = seqs[0].len();
        let mut counts = Array2::zeros((seqlen, Self::MONO_CT));
        for mut row in counts.axis_iter_mut(Axis(0)) {
            row += &pseudos;
        }

        for seq in seqs.iter() {
            if seq.len() != seqlen {
                return Err(Error::InconsistentLen);
            }

            for (mut ctr, base) in counts.axis_iter_mut(Axis(0)).zip(seq.iter()) {
                ctr += &Self::incr(*base)?;
            }
        }
        Ok(counts)
    }

    /// Returns the index of given monomer in the scores matrix using the lookup table `LK`
    /// # Arguments
    /// * `mono` - monomer, eg, b'A' for DNA or b'R' for protein
    /// # Errors
    /// * `Error::InvalidMonomer(mono)` - `mono` wasn't found in the lookup table
    fn lookup(mono: u8) -> Result<usize> {
        if mono >= 127 {
            Err(Error::InvalidMonomer { mono })
        } else {
            let idx = Self::LK[mono as usize];
            if idx == INVALID_MONO {
                Err(Error::InvalidMonomer { mono })
            } else {
                Ok(idx as usize)
            }
        }
    }

    /// Returns an array of length MONO_CT summing to 1.  used to build a PSSM
    /// from sequences, potentially including ambiguous monomers (eg, M is either A or C)
    /// # Arguments
    /// * `mono` - monomer, eg, b'A' for DNA or b'R' for protein
    /// # Errors
    /// * `Error::InvalidMonomer(mono)` - `mono` wasn't found in the lookup table
    fn incr(mono: u8) -> Result<Array1<f32>>;

    /// Returns the monomer associated with the given index; the reverse of `lookup`.
    /// Returns INVALID_MONO if the index isn't associated with a monomer.
    /// # Arguments
    /// * `idx` - the index in question
    fn rev_lk(idx: usize) -> u8;

    /// Returns the length of motif
    fn len(&self) -> usize;

    fn is_empty(&self) -> bool {
        self.len() == 0usize
    }

    /// Returns a representation of the motif using ambiguous codes.
    /// Primarily useful for DNA motifs, where ambiguous codes are
    /// common (eg, 'M' for 'A or C'); less so for proteins, where we
    /// represent any position without a dominant amino acid as an 'X'
    fn degenerate_consensus(&self) -> Vec<u8>;

    /// Accessor - returns scores matrix
    fn get_scores(&self) -> &Array2<f32>;

    /// Return sum of "worst" base at each position
    fn get_min_score(&self) -> f32;

    /// Return sum of "best" base at each position
    fn get_max_score(&self) -> f32;

    /// Returns information content of a single position.
    /// Used `info_content` method.
    /// FIXME: this should be replaced with a CTFE ... or maybe just a constant
    fn get_bits() -> f32;

    /// Returns the un-normalized sum of matching bases, useful for comparing matches from
    /// motifs of different lengths
    ///
    /// # Arguments
    /// * `seq_it` - iterator representing the query sequence
    ///
    /// # Errors
    /// * `Error::InvalidMonomer(mono)` - sequence `seq_it` contained invalid monomer `mono`
    fn raw_score<C, T>(&self, seq_it: T) -> Result<(usize, f32, Vec<f32>)>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let pssm_len = self.len();

        let mut best_start = 0;
        let mut best_score = -1.0;
        let mut best_m = Vec::new();
        // we have to look at slices, so a simple iterator won't do
        let seq = seq_it.into_iter().map(|c| *c.borrow()).collect_vec();
        let scores = self.get_scores();
        for start in 0..=seq.len() - pssm_len {
            let m: Vec<f32> = match (0..pssm_len)
                .map(|i| match Self::lookup(seq[start + i]) {
                    Err(e) => Err(e),
                    Ok(pos) => Ok(scores[[i, pos]]),
                })
                .collect()
            {
                Ok(m) => m,
                Err(e) => return Err(e),
            };
            let tot = m.iter().sum();
            if tot > best_score {
                best_score = tot;
                best_start = start;
                best_m = m;
            }
        }
        Ok((best_start, best_score, best_m))
    }

    /// Returns a `ScoredPos` struct representing the best match within the query sequence
    /// see:
    ///   MATCHTM: a tool for searching transcription factor binding sites in DNA sequences
    ///   Nucleic Acids Res. 2003 Jul 1; 31(13): 3576–3579
    ///   https://www.ncbi.nlm.nih.gov/pmc/articles/PMC169193/
    ///
    /// # Arguments
    /// * `seq_it` - iterator representing the query sequence
    ///
    /// # Errors
    /// * `Error::InvalidMonomer(mono)` - sequence `seq_it` contained invalid monomer `mono`
    /// * `Error::QueryTooShort` - sequence `seq_id` was too short
    ///
    /// # Example
    /// let pssm = DNAMotif::from_seqs(vec![
    ///            b"AAAA".to_vec(),
    ///            b"AATA".to_vec(),
    ///            b"AAGA".to_vec(),
    ///            b"AAAA".to_vec(),
    ///        ].as_ref(), None).unwrap();
    /// let start_pos = pssm.score(b"CCCCCAATA").unwrap().loc;
    fn score<C, T>(&self, seq_it: T) -> Result<ScoredPos>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let pssm_len = self.len();
        let seq = seq_it.into_iter().map(|c| *c.borrow()).collect_vec();
        if seq.len() < pssm_len {
            return Err(Error::QueryTooShort {
                motif_len: pssm_len,
                query_len: seq.len(),
            });
        }
        let min_score = self.get_min_score();
        let max_score = self.get_max_score();

        if abs_diff_eq!(max_score, min_score) {
            return Err(Error::NullMotif);
        }

        let (best_start, best_score, best_m) = self.raw_score(&seq)?;

        Ok(ScoredPos {
            loc: best_start,
            sum: (best_score - min_score) / (max_score - min_score),
            scores: best_m,
        })
    }

    /// Returns a float representing the information content of a motif; roughly the
    /// inverse of Shannon Entropy.
    /// Adapted from the information content described here:
    ///    https://en.wikipedia.org/wiki/Sequence_logo#Logo_creation
    fn info_content(&self) -> f32 {
        fn ent<'a, I>(probs: I) -> f32
        where
            I: Iterator<Item = &'a f32>,
        {
            probs
                .map(|p| {
                    if *p == 0.0 {
                        0.0
                    } else {
                        -1.0 * *p * p.log(2.0)
                    }
                })
                .sum()
        }
        let bits = Self::get_bits();
        let scores = self.get_scores();
        let mut tot = 0.0;
        for row in scores.rows() {
            tot += bits - ent(row.iter());
        }
        tot
    }
}

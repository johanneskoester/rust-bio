// Copyright 2018 Kieran Hervold
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

use super::*;
use ndarray::prelude::Array2;
use std::f32;
use std::f32::{INFINITY, NEG_INFINITY};

/// Position-specific scoring matrix for protein sequences
#[derive(Default, Clone, PartialEq, Debug)]
pub struct ProtMotif {
    /// matrix holding weights at each position, indexed by [position, base]
    pub scores: Array2<f32>,
    /// sum of "worst" base at each position
    pub min_score: f32,
    /// sum of "best" base at each position
    pub max_score: f32,
}

impl ProtMotif {
    /// Returns a Motif representing the sequences provided.
    /// # Arguments
    /// * `seqs` - sequences incorporated into motif
    /// * `pseudos` - array slice with a pseudocount for each monomer;
    ///    defaults to pssm::DEF_PSEUDO for all if None is supplied
    ///
    /// FIXME: pseudos should be an array of size MONO_CT, but that
    /// is currently impossible - see
    /// https://github.com/rust-lang/rust/issues/42863
    pub fn from_seqs(seqs: &[Vec<u8>], pseudos: Option<&[f32]>) -> Result<Self> {
        let w = Self::seqs_to_weights(seqs, pseudos)?;
        let mut m = ProtMotif {
            scores: w,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.normalize();
        m.calc_minmax();
        Ok(m)
    }

    // helper function -- normalize self.scores
    fn normalize(&mut self) {
        for i in 0..self.len() {
            let mut tot: f32 = 0.0;
            // FIXME: slices would be cleaner
            for base_i in 0..20 {
                tot += self.scores[[i, base_i]];
            }
            for base_i in 0..20 {
                self.scores[[i, base_i]] /= tot;
            }
        }
    }

    // helper function
    fn calc_minmax(&mut self) {
        let pssm_len = self.len();

        // score corresponding to sum of "worst" bases at each position
        // FIXME: iter ...
        self.min_score = 0.0;
        for i in 0..pssm_len {
            // can't use the regular min/max on f32, so we use f32::min
            let min_sc = (0..20)
                .map(|b| self.scores[[i, b]])
                .fold(INFINITY, f32::min);
            self.min_score += min_sc;
        }

        // score corresponding to "best" base at each position
        self.max_score = 0.0;
        for i in 0..pssm_len {
            let max_sc = (0..20)
                .map(|b| self.scores[[i, b]])
                .fold(NEG_INFINITY, f32::max);
            self.max_score += max_sc;
        }
    }
}

impl Motif for ProtMotif {
    const LK: [u8; 127] = [
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255,
        255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 0, 255, 4, 3, 5, 13, 7, 8, 9, 255,
        11, 10, 12, 2, 255, 14, 6, 1, 15, 16, 255, 19, 17, 255, 18, 255, 255, 255, 255, 255, 255,
        255, 0, 255, 4, 3, 5, 13, 7, 8, 9, 255, 11, 10, 12, 2, 255, 14, 6, 1, 15, 16, 255, 19, 17,
        255, 18, 255, 255, 255, 255, 255,
    ];
    const MONOS: &'static [u8] = b"ARNDCEQGHILKMFPSTWYV";
    const MONO_CT: usize = 20;

    fn rev_lk(idx: usize) -> u8 {
        if idx >= Self::MONOS.len() {
            INVALID_MONO
        } else {
            Self::MONOS[idx]
        }
    }
    fn incr(mono: u8) -> Result<Array1<f32>> {
        if mono >= 127 {
            Err(Error::InvalidMonomer { mono })
        } else {
            if mono == b'X' {
                Ok(Array1::from_elem(Self::MONO_CT, 1.0 / Self::MONO_CT as f32))
            } else {
                let idx = Self::LK[mono as usize] as usize;
                let mut v = Array1::zeros(Self::MONO_CT);
                v[idx] += 1.0;
                Ok(v)
            }
        }
    }
    fn len(&self) -> usize {
        self.scores.dim().0
    }
    fn get_scores(&self) -> &Array2<f32> {
        &self.scores
    }
    fn get_min_score(&self) -> f32 {
        self.min_score
    }
    fn get_max_score(&self) -> f32 {
        self.max_score
    }
    fn get_bits() -> f32 {
        20f32.log2()
    }
    fn degenerate_consensus(&self) -> Vec<u8> {
        let len = self.len();
        let mut res = Vec::with_capacity(len);
        for pos in 0..len {
            let mut fracs = (0..20)
                .map(|b| (self.scores[[pos, b]], b))
                .collect::<Vec<(f32, usize)>>();
            // note: reverse sort
            fracs.sort_by(|a, b| b.partial_cmp(a).unwrap());

            res.push(if fracs[0].0 > 0.5 && fracs[0].0 > 2.0 * fracs[1].0 {
                Self::MONOS[fracs[0].1]
            } else {
                b'X'
            });
        }
        res
    }
}

/// Return a ProtMotif wrapping an Array2 representing amino acid
/// weights at each position.  The dimensions and contents of this
/// array are unchecked, and it is incumbent on the user to ensure
/// the correct dimensions are used (ie, SEQ_LEN x 20), and no zeros
/// appear in the array.
impl From<Array2<f32>> for ProtMotif {
    fn from(scores: Array2<f32>) -> Self {
        let mut m = ProtMotif {
            scores,
            min_score: 0.0,
            max_score: 0.0,
        };
        m.normalize();
        m.calc_minmax();
        m
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array;

    #[test]
    fn test_info_content() {
        let pssm = ProtMotif::from_seqs(vec![b"AAAA".to_vec()].as_ref(), Some(&[0.0; 20])).unwrap();
        assert_relative_eq!(
            pssm.info_content(),
            ProtMotif::get_bits() * 4.0,
            epsilon = f32::EPSILON
        );
    }

    #[test]
    fn test_scoring() {
        // should match "ARND"
        let m: Array2<f32> = Array::from(vec![
            0.81, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.81, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            0.81, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.81, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
            0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
        ])
        .into_shape((4, 20))
        .unwrap();
        let pssm: ProtMotif = m.into();
        let scored_pos = pssm.score(b"AAAAARNDAAA").unwrap();
        assert_eq!(scored_pos.loc, 4);
    }

    #[test]
    fn test_mono_err() {
        let pssm = ProtMotif::from_seqs(&[b"ARGN".to_vec()], None).unwrap();
        assert!(matches!(
            pssm.score(b"AAAABAAAAAAAAA"),
            Err(Error::InvalidMonomer { mono: b'B' })
        ));
    }

    #[test]
    fn test_inconsist_err() {
        assert!(matches!(
            ProtMotif::from_seqs(
                &[b"NNNNN".to_vec(), b"RRRRR".to_vec(), b"C".to_vec()],
                Some(&[0.0; 20])
            ),
            Err(Error::InconsistentLen)
        ));
    }

    #[test]
    fn test_degenerate_consensus_same_bases() {
        let pssm = ProtMotif::from_seqs(
            &[b"QVTYNDSA".to_vec(), b"QVTYNDSA".to_vec()],
            Some(&[0.0; 20]),
        )
        .unwrap();
        assert_eq!(pssm.degenerate_consensus(), b"QVTYNDSA".to_vec());
    }

    #[test]
    fn test_degenerate_consensus_x() {
        let pssm = ProtMotif::from_seqs(
            &[b"QVTYNDSA".to_vec(), b"ASDNYTVQ".to_vec()],
            Some(&[0.0; 20]),
        )
        .unwrap();
        assert_eq!(pssm.degenerate_consensus(), b"XXXXXXXX".to_vec());
    }
    fn test_degenerate_input() {
        let pssm = ProtMotif::from_seqs(
            &[b"AAAAARNDAAA".to_vec(), b"AAAAARNDXAA".to_vec()],
            Some(&[0.0; 20]),
        )
        .unwrap();
        assert_eq!(pssm.degenerate_consensus(), b"AAAAARNDXAA".to_vec());
    }
}

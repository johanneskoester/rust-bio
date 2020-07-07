// Copyright 2014-2015 Johannes Köster, Peer Aramillo Irizar.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Implementation of alphabets and useful utilities.
//!
//! # Example
//!
//! ```rust
//! use bio::alphabets;
//! let alphabet = alphabets::dna::alphabet();
//! assert!(alphabet.is_word(b"AACCTgga"));
//! assert!(!alphabet.is_word(b"AXYZ"));
//! ```

use std::borrow::Borrow;
use std::mem;

use bit_set::BitSet;
use vec_map::VecMap;

pub mod dna;
pub mod protein;
pub mod rna;

pub type SymbolRanks = VecMap<u8>;

/// Representation of an alphabet.
pub struct Alphabet {
    pub symbols: BitSet,
}

impl Alphabet {
    /// Create new alphabet from given symbols.
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// assert!(dna_alphabet.is_word(b"GAttACA"));
    /// ```
    pub fn new<C, T>(symbols: T) -> Self
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let mut s = BitSet::new();
        s.extend(symbols.into_iter().map(|c| *c.borrow() as usize));

        Alphabet { symbols: s }
    }

    /// Insert symbol into alphabet.
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let mut dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// assert!(!dna_alphabet.is_word(b"N"));
    /// dna_alphabet.insert(78);
    /// assert!(dna_alphabet.is_word(b"N"));
    /// ```
    pub fn insert(&mut self, a: u8) {
        self.symbols.insert(a as usize);
    }

    /// Check if given text is a word over the alphabet.
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// assert!(dna_alphabet.is_word(b"GAttACA"));
    /// assert!(!dna_alphabet.is_word(b"42"));
    /// ```
    pub fn is_word<C, T>(&self, text: T) -> bool
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .all(|c| self.symbols.contains(*c.borrow() as usize))
    }

    /// Return lexicographically maximal symbol.
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// assert_eq!(dna_alphabet.max_symbol(), Some(116));  // max symbol is "t"
    /// let empty_alphabet = alphabets::Alphabet::new(b"");
    /// assert_eq!(empty_alphabet.max_symbol(), None);
    /// ```
    pub fn max_symbol(&self) -> Option<u8> {
        self.symbols.iter().max().map(|a| a as u8)
    }

    /// Return size of the alphabet.
    ///
    /// Upper and lower case representations of the same character
    /// are counted as distinct characters.
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// assert_eq!(dna_alphabet.len(), 8);
    /// ```
    pub fn len(&self) -> usize {
        self.symbols.len()
    }

    /// Is this alphabet empty?
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// assert!(!dna_alphabet.is_empty());
    /// let empty_alphabet = alphabets::Alphabet::new(b"");
    /// assert!(empty_alphabet.is_empty());
    /// ```
    pub fn is_empty(&self) -> bool {
        self.symbols.is_empty()
    }
}

/// Tools based on transforming the alphabet symbols to their lexicographical ranks.
///
/// Lexicographical rank is computed using `u8` representations,
/// i.e. ASCII codes, of the input characters.
#[derive(Serialize, Deserialize)]
pub struct RankTransform {
    pub ranks: SymbolRanks,
}

impl RankTransform {
    /// Construct a new `RankTransform`.
    pub fn new(alphabet: &Alphabet) -> Self {
        let mut ranks = VecMap::new();
        for (r, c) in alphabet.symbols.iter().enumerate() {
            ranks.insert(c, r as u8);
        }

        RankTransform { ranks }
    }

    /// Get the rank of symbol `a`.
    pub fn get(&self, a: u8) -> u8 {
        *self.ranks.get(a as usize).expect("Unexpected character.")
    }

    /// Transform a given `text`.
    pub fn transform<C, T>(&self, text: T) -> Vec<u8>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        text.into_iter()
            .map(|c| {
                *self
                    .ranks
                    .get(*c.borrow() as usize)
                    .expect("Unexpected character in text.")
            })
            .collect()
    }

    /// Iterate over q-grams (substrings of length q) of given `text`. The q-grams are encoded
    /// as `usize` by storing the symbol ranks in log2(|A|) bits (with |A| being the alphabet size).
    ///
    /// If q is larger than usize::BITS / log2(|A|), this method fails with an assertion.
    pub fn qgrams<C, T>(&self, q: u32, text: T) -> QGrams<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!(
            (bits * q) as usize <= mem::size_of::<usize>() * 8,
            "Expecting q to be smaller than usize / log2(|A|)"
        );

        let mut qgrams = QGrams {
            text: text.into_iter(),
            ranks: self,
            bits,
            mask: (1 << (q * bits)) - 1,
            qgram: 0,
        };

        for _ in 0..q - 1 {
            qgrams.next();
        }

        qgrams
    }

    /// Restore alphabet from transform.
    pub fn alphabet(&self) -> Alphabet {
        let mut symbols = BitSet::with_capacity(self.ranks.len());
        symbols.extend(self.ranks.keys());
        Alphabet { symbols }
    }
}

/// Iterator over q-grams.
pub struct QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    text: T,
    ranks: &'a RankTransform,
    bits: u32,
    mask: usize,
    qgram: usize,
}

impl<'a, C, T> QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    /// Push a new character into the current qgram.
    fn qgram_push(&mut self, a: u8) {
        self.qgram <<= self.bits;
        self.qgram |= a as usize;
        self.qgram &= self.mask;
    }
}

impl<'a, C, T> Iterator for QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        match self.text.next() {
            Some(a) => {
                let b = self.ranks.get(*a.borrow());
                self.qgram_push(b);
                Some(self.qgram)
            }
            None => None,
        }
    }
}

#[cfg(tests)]
mod tests {
    #[test]
    fn test_serde() {
        use serde::{Deserialize, Serialize};
        fn impls_serde_traits<S: Serialize + Deserialize>() {}

        impls_serde_traits::<RankTransform>();
    }
}

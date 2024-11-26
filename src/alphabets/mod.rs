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

use bit_set::BitSet;
use vec_map::VecMap;

pub mod dna;
pub mod protein;
pub mod rna;

pub type SymbolRanks = VecMap<u8>;

/// Representation of an alphabet.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug)]
pub struct Alphabet {
    pub symbols: BitSet,
}

impl Alphabet {
    /// Create new alphabet from given symbols.
    ///
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// // Create an alphabet (note that a DNA alphabet is already available in bio::alphabets::dna).
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// // Check whether a given text is a word over the alphabet.
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
    /// Complexity: O(1)
    ///
    /// # Example
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
    /// Complexity: O(n), where n is the length of the text.
    ///
    /// # Example
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
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// assert_eq!(dna_alphabet.max_symbol(), Some(116)); // max symbol is "t"
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
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
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
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
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

    /// Return a new alphabet taking the intersect between this and other.
    ///
    /// # Example
    /// ```
    /// use bio::alphabets;
    ///
    /// let alpha_a = alphabets::Alphabet::new(b"acgtACGT");
    /// let alpha_b = alphabets::Alphabet::new(b"atcgMVP");
    /// let intersect_alpha = alpha_a.intersection(&alpha_b);
    ///
    /// assert_eq!(intersect_alpha, alphabets::Alphabet::new(b"atcg"));
    /// ```
    pub fn intersection(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.intersection(&other.symbols).collect(),
        };
    }

    /// Return a new alphabet taking the difference between this and other.
    ///
    /// # Example
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// let dna_alphabet_upper = alphabets::Alphabet::new(b"ACGT");
    /// let dna_lower = dna_alphabet.difference(&dna_alphabet_upper);
    ///
    /// assert_eq!(dna_lower, alphabets::Alphabet::new(b"atcg"));
    /// ```
    pub fn difference(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.difference(&other.symbols).collect(),
        };
    }

    /// Return a new alphabet taking the union between this and other.
    ///
    /// # Example
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ATCG");
    /// let tokenize_alpha = alphabets::Alphabet::new(b"?|");
    /// let alpha = dna_alphabet.union(&tokenize_alpha);
    ///
    /// assert_eq!(alpha, alphabets::Alphabet::new(b"ATCG?|"));
    /// ```
    pub fn union(&self, other: &Alphabet) -> Self {
        return Alphabet {
            symbols: self.symbols.union(&other.symbols).collect(),
        };
    }
}

/// Tools based on transforming the alphabet symbols to their lexicographical ranks.
///
/// Lexicographical rank is computed using `u8` representations,
/// i.e. ASCII codes, of the input characters.
/// For example, assuming that the alphabet consists of the symbols `A`, `C`, `G`, and `T`, this
/// will yield ranks `0`, `1`, `2`, `3` for them, respectively.
///
/// `RankTransform` can be used in to perform bit encoding for texts over a
/// given alphabet via `bio::data_structures::bitenc`.
#[derive(Default, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize, Deserialize)]
pub struct RankTransform {
    pub ranks: SymbolRanks,
}

impl RankTransform {
    /// Construct a new `RankTransform`.
    ///
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    /// ```
    pub fn new(alphabet: &Alphabet) -> Self {
        let mut ranks = VecMap::new();
        for (r, c) in alphabet.symbols.iter().enumerate() {
            ranks.insert(c, r as u8);
        }

        RankTransform { ranks }
    }

    /// Get the rank of symbol `a`.
    ///
    /// This method panics for characters not contained in the alphabet.
    ///
    /// Complexity: O(1)
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    /// assert_eq!(dna_ranks.get(65), 0); // "A"
    /// assert_eq!(dna_ranks.get(116), 7); // "t"
    /// ```
    #[inline]
    pub fn get(&self, a: u8) -> u8 {
        *self.ranks.get(a as usize).expect("Unexpected character.")
    }

    /// Transform a given `text` into a vector of rank values.
    ///
    /// Complexity: O(n), where n is the length of the text.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    /// let text = b"aAcCgGtT";
    /// assert_eq!(dna_ranks.transform(text), vec![4, 0, 5, 1, 6, 2, 7, 3]);
    /// ```
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
    ///
    /// Complexity: O(n), where n is the length of the text.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    ///
    /// let q_grams: Vec<usize> = dna_ranks.qgrams(2, b"ACGT").collect();
    /// assert_eq!(q_grams, vec![1, 10, 19]);
    /// ```
    pub fn qgrams<C, T>(&self, q: u32, text: T) -> QGrams<'_, C, T::IntoIter>
    where
        C: Borrow<u8>,
        T: IntoIterator<Item = C>,
    {
        assert!(q > 0, "Expecting q-gram length q to be larger than 0.");
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!(
            bits * q <= usize::BITS,
            "Expecting q to be smaller than usize / log2(|A|)"
        );

        let mut qgrams = QGrams {
            text: text.into_iter(),
            ranks: self,
            q,
            bits,
            mask: 1usize.checked_shl(q * bits).unwrap_or(0).wrapping_sub(1),
            qgram: 0,
        };

        for _ in 0..q - 1 {
            qgrams.next();
        }

        qgrams
    }

    /// Iterate over q-grams (substrings of length q) of given `text` in reverse. The q-grams are encoded
    /// as `usize` by storing the symbol ranks in log2(|A|) bits (with |A| being the alphabet size).
    ///
    /// If q is larger than usize::BITS / log2(|A|), this method fails with an assertion.
    ///
    /// Complexity: O(n), where n is the length of the text.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGTacgt");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    ///
    /// let q_grams: Vec<usize> = dna_ranks.rev_qgrams(2, b"ACGT").collect();
    /// assert_eq!(q_grams, vec![19, 10, 1]);
    /// ```
    pub fn rev_qgrams<C, IT, T>(&self, q: u32, text: IT) -> RevQGrams<'_, C, T>
    where
        C: Borrow<u8>,
        T: DoubleEndedIterator<Item = C>,
        IT: IntoIterator<IntoIter = T>,
    {
        assert!(q > 0, "Expecting q-gram length q to be larger than 0.");
        let bits = (self.ranks.len() as f32).log2().ceil() as u32;
        assert!(
            bits * q <= usize::BITS,
            "Expecting q to be smaller than usize / log2(|A|)"
        );

        let mut rev_qgrams = RevQGrams {
            text: text.into_iter(),
            ranks: self,
            q,
            bits,
            left_shift: (q - 1) * bits,
            qgram: 0,
        };

        for _ in 0..q - 1 {
            rev_qgrams.next();
        }

        rev_qgrams
    }

    /// Restore alphabet from transform.
    ///
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"acgtACGT");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    /// assert_eq!(dna_ranks.alphabet().symbols, dna_alphabet.symbols);
    /// ```
    pub fn alphabet(&self) -> Alphabet {
        let mut symbols = BitSet::with_capacity(self.ranks.len());
        symbols.extend(self.ranks.keys());
        Alphabet { symbols }
    }

    /// Compute the number of bits required to encode the largest rank value.
    ///
    /// For example, the alphabet `b"ACGT"` with 4 symbols has the maximal rank
    /// 3, which can be encoded in 2 bits.
    ///
    /// This value can be used to create a `data_structures::bitenc::BitEnc`
    /// bit encoding tailored to the given alphabet.
    ///
    /// Complexity: O(n), where n is the number of symbols in the alphabet.
    ///
    /// # Example
    ///
    /// ```
    /// use bio::alphabets;
    ///
    /// let dna_alphabet = alphabets::Alphabet::new(b"ACGT");
    /// let dna_ranks = alphabets::RankTransform::new(&dna_alphabet);
    /// assert_eq!(dna_ranks.get_width(), 2);
    /// let dna_n_alphabet = alphabets::Alphabet::new(b"ACGTN");
    /// let dna_n_ranks = alphabets::RankTransform::new(&dna_n_alphabet);
    /// assert_eq!(dna_n_ranks.get_width(), 3);
    /// ```
    pub fn get_width(&self) -> usize {
        (self.ranks.len() as f32).log2().ceil() as usize
    }
}

/// Iterator over q-grams.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: Iterator<Item = C>,
{
    text: T,
    ranks: &'a RankTransform,
    q: u32,
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
    #[inline]
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

    #[inline]
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

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.text.size_hint()
    }
}

impl<'a, C, T> ExactSizeIterator for QGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: ExactSizeIterator<Item = C>,
{
}

/// Iterator over q-grams.
#[derive(Copy, Clone, Eq, PartialEq, Ord, PartialOrd, Hash, Debug, Serialize)]
pub struct RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C>,
{
    text: T,
    ranks: &'a RankTransform,
    q: u32,
    bits: u32,
    left_shift: u32,
    qgram: usize,
}

impl<'a, C, T> RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C>,
{
    /// Push a new character into the current qgram in the reverse direction.
    #[inline]
    fn qgram_push_rev(&mut self, a: u8) {
        self.qgram >>= self.bits;
        self.qgram |= (a as usize) << self.left_shift;
    }
}

impl<'a, C, T> Iterator for RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C>,
{
    type Item = usize;

    #[inline]
    fn next(&mut self) -> Option<usize> {
        match self.text.next_back() {
            Some(a) => {
                let b = self.ranks.get(*a.borrow());
                self.qgram_push_rev(b);
                Some(self.qgram)
            }
            None => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.text.size_hint()
    }
}

impl<'a, C, T> ExactSizeIterator for RevQGrams<'a, C, T>
where
    C: Borrow<u8>,
    T: DoubleEndedIterator<Item = C> + ExactSizeIterator<Item = C>,
{
}

/// Returns the english ascii lower case alphabet.
pub fn english_ascii_lower_alphabet() -> Alphabet {
    Alphabet::new(&b"abcdefghijklmnopqrstuvwxyz"[..])
}

/// Returns the english ascii upper case alphabet.
pub fn english_ascii_upper_alphabet() -> Alphabet {
    Alphabet::new(&b"ABCDEFGHIJKLMNOPQRSTUVWXYZ"[..])
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_alphabet_eq() {
        assert_eq!(Alphabet::new(b"ATCG"), Alphabet::new(b"ATCG"));
        assert_eq!(Alphabet::new(b"ATCG"), Alphabet::new(b"TAGC"));
        assert_ne!(Alphabet::new(b"ATCG"), Alphabet::new(b"ATC"));
    }

    /// When `q * bits == usize::BITS`, make sure that `1<<(q*bits)` does not overflow.
    #[test]
    fn test_qgram_shiftleft_overflow() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTG".repeat(100);
        transform.qgrams(usize::BITS / 2, text);
    }

    #[test]
    fn test_qgrams() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTGA";
        let mut qgrams = transform.qgrams(4, text);
        assert_eq!(qgrams.next(), transform.qgrams(4, b"ACTG").next());
        assert_eq!(qgrams.next(), transform.qgrams(4, b"CTGA").next());
        assert_eq!(qgrams.next(), None);
    }

    #[test]
    fn test_rev_qgrams() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTGA";
        let mut rev_qgrams = transform.rev_qgrams(4, text);
        assert_eq!(rev_qgrams.next(), transform.qgrams(4, b"CTGA").next());
        assert_eq!(rev_qgrams.next(), transform.qgrams(4, b"ACTG").next());
        assert_eq!(rev_qgrams.next(), None);
    }

    #[test]
    fn test_exactsize_iterator() {
        let alphabet = Alphabet::new(b"ACTG");
        let transform = RankTransform::new(&alphabet);
        let text = b"ACTGACTG";

        let mut qgrams = transform.qgrams(4, text);
        assert_eq!(qgrams.len(), 5);
        qgrams.next();
        assert_eq!(qgrams.len(), 4);

        let mut rev_qgrams = transform.rev_qgrams(4, text);
        assert_eq!(rev_qgrams.len(), 5);
        rev_qgrams.next();
        assert_eq!(rev_qgrams.len(), 4);

        let qgrams = transform.qgrams(4, b"AC");
        assert_eq!(qgrams.len(), 0);
        let rev_qgrams = transform.rev_qgrams(4, b"AC");
        assert_eq!(rev_qgrams.len(), 0);
    }
}

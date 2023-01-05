// Copyright 2018 Johannes Köster.
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! Fast lookup table for arbitrary floating point functions.
//! # Examples
//! ## Easy:
//! ```
//! use bio::data_structures::interpolation_table::*;
//!
//! let table = InterpolationTable::new(0.0, 10.0, 5, |x| x.powf(2.0));
//! assert_eq!(table.get(3.0), 9.0);
//! assert_eq!(table.get(5.0), 25.0);
//! ```
//! ## More complicated:
//! ```
//! extern crate approx;
//! fn main() {
//!     use bio::data_structures::interpolation_table::*;
//!     use approx::assert_relative_eq;
//!
//!     let table = InterpolationTable::new(0.0, 10.0, 5, |x| x.ln_1p());
//!     for &x in &[0.02, 0.04, 0.45678686, 0.23875, 1.45345e-6] {
//!         assert_relative_eq!(table.get(x), x.ln_1p(), epsilon = 0.00001);
//!     }
//! }

#[inline]
pub fn interpolate(a: f64, b: f64, fraction: f64) -> f64 {
    a * (1.0 - fraction) + b * fraction
}

/// Fast lookup table for arbitrary floating point functions.
/// This can be used to e.g., provide fast lookups of distribution values.
/// Input values are sampled with a given precision and results are stored in a vector.
/// During lookup, infimum and supremum of a given value are calculated and the result is
/// interpolated.
#[derive(Default, Clone, PartialEq, PartialOrd, Debug, Serialize, Deserialize)]
pub struct InterpolationTable<F: Fn(f64) -> f64> {
    inner: Vec<f64>,
    func: F,
    offset: usize,
    min_x: f64,
    max_x: f64,
    shift: f64,
}

impl<F: Fn(f64) -> f64> InterpolationTable<F> {
    /// Create a new `InterpolationTable`.
    ///
    /// # Arguments
    ///
    /// * `min_x` - minimum sample value
    /// * `max_x` - maximum sample value
    /// * `frac_digits` - number of fraction digits to store in sample
    /// * `func` - Function to emulate.
    ///
    /// If given value is outside of min_x and max_x, the lookup falls back to applying the
    /// function itself.
    /// The table size grows with the number of fraction digits.
    /// Space Complexity: O(m * 10^n), where `m = max_x - min_x` and `n = frac_digits`
    pub fn new(min_x: f64, max_x: f64, frac_digits: i32, func: F) -> Self {
        let shift = 10.0_f64.powi(frac_digits);
        let offset = (min_x * shift) as usize;
        let mut table = InterpolationTable {
            inner: Vec::new(),
            func,
            min_x,
            max_x,
            shift,
            offset,
        };

        for i in table.index(min_x)..table.index(max_x) {
            let x = i as f64 / shift;
            table.inner.push((table.func)(x));
        }

        table
    }

    /// Return lower bound index for given f64.
    #[inline]
    fn index(&self, x: f64) -> usize {
        (x * self.shift) as usize - self.offset
    }

    /// Lookup given value in table, and interpolate the result between the sampled values if
    /// necessary. This provides an approximation that is better the more fraction digits are
    /// used to generate this table.
    /// Time Complexity for lookup: O(1) if `min_x <= x < max_x` and O(func(x)) otherwise.
    pub fn get(&self, x: f64) -> f64 {
        if x < self.min_x || x >= self.max_x {
            (self.func)(x)
        } else {
            let i = self.index(x);
            // interpolate
            let fraction = (x * self.shift - i as f64) / self.shift;
            interpolate(self.inner[i], self.inner[i + 1], fraction)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interpolation_table() {
        let table = InterpolationTable::new(0.0, 10.0, 5, |x| x.ln_1p());

        for &x in &[0.02, 0.04, 0.45678686, 0.23875, 1.45345e-6] {
            assert_relative_eq!(table.get(x), x.ln_1p(), epsilon = 0.00001);
        }
    }
}

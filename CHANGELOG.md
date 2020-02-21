# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

# [0.30.0] - 2019-11-14
- Bayesian models now allow to access internals.
- Various small bug fixes.

# [0.29.0] - 2019-09-27
- Migrate error handling to the snafu crate (this is an API breaking change).
- Fix edge cases in pairwise alignment.
- Fix error in backward search if symbol isn't found.

# [0.28.1] - 2019-06-28
- Fix select in RankSelect in cases where many superblocks have the same rank.

# [0.28.0] - 2019-06-19
- Myers bit-parallel pattern matching now supports arbitrarily long patterns via bit vectors (thanks to @markschl).
- Minor documentation updates (thanks to @anders-was-here).

# [0.27.0] - 2019-05-31
- Implement sequence-read-trait for FASTQ records.
- Cleanup dependencies.

# [0.26.1] - 2019-05-10
- Fix a bug in `select_1` and `select_0` that would lead to too large answers.

# [0.26.0] - 2019-05-09
- Added a trait system for computing Bayesian statistical models.
- Added an implementation of MSA via partial order alignment.
- Performance improvements to FASTQ reader.

# [0.25.0] - 2018-12-12
- Added `FQRead` and `FARead` traits to `FastaReader` and `FastqReader` to be more flexible with input types. This allows to use readers on gzipped and on plain text input interchangeably.
- Added an implementation of Bayes Factors and evidence scoring using the method of Kass and Raftery.

# [0.24.0] - 2018-11-26
- API overhaul to become more flexible when accepting text iterators. Now, anything that iterates over something can be borrowed as u8 is allowed.
- FMIndex and FMDIndex now also allow plain owned versions of BWT, Less and Occ. This should greatly simplify their usage.
- PairHMM and LogProb implementation has seen extensive performance improvements. Among that, (a) the usage of a fast approximation of exp() as presented by [Kopczynsi 2017](https://eldorado.tu-dortmund.de/bitstream/2003/36203/1/Dissertation_Kopczynski.pdf), and (b) banding of the pairHMM matrix with a given maximum edit distance.
- All IO records now support serde.

# [0.23.0] - 2018-11-06
- Generalized Myers pattern matching algorithm to arbitrary unsigned integer types (u64, u128) (thanks to @markschl).
- Implemented optional traceback and alignment output for Myers pattern matching algorithm (thanks to @markschl).
- Use Strand type from bio-types crate in BED module (thanks to @ingolia).
- Added an IntervalTree based data structure for looking up overlaps between annotation types (thanks to @ingolia).
- Various bug fixes.

# [0.22.0] - 2018-08-01
- Added HMM implementation (thanks to @holtgrewe).
- Moved Alignment types to `bio_types` crate (thanks to @pmarks).
- Ignore comment lines in GTF/GFF files (thanks to Yasunobu Okamura).
- API usability improvements.

# [0.21.0] - 2018-06-19
- Added PSSM implementation (thanks to @hervold).

# [0.20.0] - 2018-06-01
- Refactored RankSelect API to consistently use u64.
- Use bv crate in suffix array implementation.

# [0.19.0] - 2018-05-25
- rank-0 and select-0 in RankSelect.
- use bv crate for RankSelect.

## [0.18.0] - 2018-05-04
- More flexible FASTA API.
- Fixed bug in KMP.

## [0.17.0] - 2018-02-22
- Bug fix in Ukkonen algorithm
- Convenience improvements to API

## [0.16.0] - 2018-01-05
- Pairwise alignment has been rewritten to support banded alignment and clips.
- Various minor API additions and improvements.
- Several small bug fixes. 

## [0.15.0] - 2017-11-20
- Add pair hidden markov model implementation to calculate the probability of two sequences being related.
- Various minor bug fixes and usability improvements.

## [0.14.2] - 2017-08-30
- Improved numerical stability of CDF construction.
- Speed improvements to occurrence array lookups in FM-index.
- Improved GFF/GTF variant format handling.
- Improved robustness of credible interval calculating in CDF.
- Bug fixes for log probability implementation.

## [0.14.1] - 2017-06-23
### Changed
- Replace nalgebra dependency with ndarray crate.

## [0.14.0] - 2017-06-15
### Changed
- GTF/GFF reader can now handle duplicate keys.
- Updated dependencies.
- RNA alphabet.
- Improved FASTQ reader.
- Fixes in alignment algorithm.


## [0.13.0] - 2017-05-09
### Changed
- fasta::IndexedReader now also provides an iterator.
- IntervalTree provides a mutable iterator.
- Various fixes to Fasta IO.
- Fixed calculation of expected FDR.

## [0.12.0] - 2017-04-03
### Changed
- Improved distance API.
- Moved Strand into utils.
- More robust gff/gtf parsing.


## [0.11.0] - 2017-02-16
### Changed
- Improved IntervalTree API.
- Updated dependencies.
- Speed improvements in alignment module.
- Improved test coverage.
- Speed improvements in fmindex module.

## [0.10.0] - 2016-11-02
### Added
- An interval tree implementation.
- Initial utilities for bayesian statistics.
### Changed
- Various small improvements to log-space probability API.

## [0.9.0] - 2016-08-18
### Added
- Implementation of discrete probability distributions via cumulative distribution functions.
### Changed
- Log-space probabilities have been refactored into newtypes.
- Performance improvements for FMIndex implementation.
- Improved documentation.

## [0.8.0] - 2016-07-20
### Changed
- Writers in the io module no longer take ownership of the given record.
- Various cosmetic changes.

## [0.7.0] - 2016-07-06
### Changed
- Reverse complement API has been refactored into plain functions.
- Reverse complement now supports the whole IUPAC alphabet.
- Various algorithms take now IntoTextIterator instead of only slices.
- Fasta reader and writer treat sequence names as strings.
- Refactoring of suffix array + fmindex API to provide more flexibility.

## [0.6.0] - 2016-05-09
### Changed
- Type aliases for various text representations.
- Pattern matching algorithms take both iterators and slices where possible.
- logprobs::cumsum has been refactored to return an iterator.
- support for subtraction of logprobs.

## [0.5.0] - 2016-02-24
### Added
- Support for [serde](https://github.com/serde-rs/serde) serialization when used in combination with rust nightly (@dikaiosune).

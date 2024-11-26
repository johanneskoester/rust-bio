//! Readers and writers for common bioinformatics file formats.

pub mod bed;
pub mod fasta;
pub mod fastq;
pub mod gaf;
pub mod gff;
#[cfg(feature = "phylogeny")]
pub mod newick;

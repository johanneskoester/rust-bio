#![feature(test)]

extern crate test;

use test::Bencher;

use bio::io::fasta;

use std::io::Seek;
use std::io::Write;

use rand::SeedableRng;

const NUCS: [u8; 4] = [b'A', b'C', b'T', b'G'];

fn random_seq<RNG>(length: usize, rng: &mut RNG) -> Vec<u8>
where
    RNG: rand::Rng,
{
    (0..length).map(|_| NUCS[rng.gen_range(0..=3)]).collect()
}

fn write_sequence(file: &mut std::fs::File, seed: u64) {
    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);

    for i in 0..1000 {
        writeln!(
            file,
            ">{}\n{}",
            i,
            String::from_utf8(random_seq(300, &mut rng)).unwrap()
        )
        .unwrap();
    }
}

mod fasta_buffer {
    use std::io;

    use super::*;

    #[bench]
    fn default(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            fasta::Reader::new(&mut tempfile).records().for_each(|_| ())
        });
    }

    #[bench]
    fn wrapped_default(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            fasta::Reader::new(io::BufReader::new(&mut tempfile))
                .records()
                .for_each(|_| ())
        });
    }

    #[bench]
    fn capacity_default(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            fasta::Reader::with_capacity(8192, &mut tempfile)
                .records()
                .for_each(|_| ())
        });
    }

    #[bench]
    fn bufread_default(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            let buffer = io::BufReader::new(&mut tempfile);
            fasta::Reader::from_bufread(buffer)
                .records()
                .for_each(|_| ())
        });
    }

    #[bench]
    fn wrapped_32768(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            fasta::Reader::new(io::BufReader::with_capacity(32768, &mut tempfile))
                .records()
                .for_each(|_| ())
        });
    }

    #[bench]
    fn capacity_32768(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            fasta::Reader::with_capacity(32768, &mut tempfile)
                .records()
                .for_each(|_| ())
        });
    }

    #[bench]
    fn bufread_32768(b: &mut Bencher) {
        let mut tempfile = tempfile::tempfile().unwrap();

        write_sequence(&mut tempfile, 42);

        b.iter(|| {
            tempfile.seek(io::SeekFrom::Start(0)).unwrap();
            let buffer = io::BufReader::with_capacity(32768, &mut tempfile);
            fasta::Reader::from_bufread(buffer)
                .records()
                .for_each(|_| ())
        });
    }
}

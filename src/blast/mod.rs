
//! Blastn search from NCBI with data formating related
//! 
//! Transformation and transcription from dna and rna strand, between other formating data functions
//! ```
//! use bio::blast;
//! let dna_string = vec![b'C',b'T',b'A'];
//! let rna_string = blast::transcription(dna_string);
//! assert!(rna_string == vec![b'C',b'U',b'A']);
  
//! let prot_string = blast::traduction(rna_string);
//! assert!(prot_string == vec![b'L']);
//! ```

use crate::utils::Text;
//use vec;

pub use self::blast_req::blastn_req;

pub mod blast_req;


/// Transcription of a dna (Text :type) into the rna version
pub fn transcription(dna_str : Text) -> Text {
    //let alphab = alphabets::dna::alphabet();
    let mut rna_str : Text = vec![];
    for i in &dna_str {
        let nuc_dna = i;
        match nuc_dna{
            84=>rna_str.push(85),
            116=>rna_str.push(117),
            _=>rna_str.push(*nuc_dna),
        }
    }
    rna_str
}

/// traduction of a rna (Text :type) into the aminoacid strand
pub fn traduction(rna_str : Text) -> Text{
    let mut codo : Text = vec![];
    let mut protein_str : Text = vec![];
    for i in &rna_str {
        if codo.len() < 3 {
            codo.push(*i);
        }
        if codo.len() == 3{
            match codo[0]{
                b'C'=> match codo[1]{
                    b'A'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'H'),
                        _=>protein_str.push(b'Q'),
                    },
                    b'U'=>protein_str.push(b'L'),
                    b'C'=>protein_str.push(b'P'),
                    _=>protein_str.push(b'R'), //G
                },
                b'G'=> match codo[1]{
                    b'A'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'D'),
                        _=>protein_str.push(b'E'),
                    },
                    b'U'=>protein_str.push(b'V'),
                    b'C'=>protein_str.push(b'A'),
                    _=>protein_str.push(b'G'), //G
                },
                b'A'=> match codo[1]{
                    b'A'=> match codo[2]{
                        b'G' | b'A' =>protein_str.push(b'K'),
                        _=>protein_str.push(b'N'),
                    },
                    b'G'=> match codo[2]{
                        b'G' | b'A' =>protein_str.push(b'R'),
                        _=>protein_str.push(b'S'),
                    },
                    b'U'=>match codo[2]{
                        b'G' =>protein_str.push(b'M'),
                        _=>protein_str.push(b'I'),
                    },
                    _=>protein_str.push(b'P'), //C
                },
                b'U'=> match codo[1]{
                    b'U'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'F'),
                        _=>protein_str.push(b'L'),
                    },
                    b'A'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'Y'),
                        _=>protein_str.push(b'O'),
                    },
                    b'G'=> match codo[2]{
                        b'U' | b'C' =>protein_str.push(b'C'),
                        b'G'=>protein_str.push(b'W'),
                        _=>protein_str.push(b'O'),
                    },
                    _=>protein_str.push(b'S'), //C
                },
                _=>protein_str.push(b'X'),
            }

            codo.clear();
        }
    }
    
    protein_str
}

///Transformation from uniletter format to three letter format of a protein strand in Text vector
pub fn protein_to_three(protein_str : Text) -> Vec<String>{
    let mut three_prot : Vec<String> = vec![];
    for i in protein_str {
        match i{
            b'F'=>three_prot.push("Phe".to_string()),
            b'L'=>three_prot.push("Leu".to_string()),
            b'S'=>three_prot.push("Ser".to_string()),
            b'Y'=>three_prot.push("Tyr".to_string()),
            b'O'=>three_prot.push("Stop".to_string()),
            b'C'=>three_prot.push("Cys".to_string()),
            b'W'=>three_prot.push("Trp".to_string()),
            b'P'=>three_prot.push("Pro".to_string()),
            b'H'=>three_prot.push("His".to_string()),
            b'Q'=>three_prot.push("Gln".to_string()),
            b'R'=>three_prot.push("Arg".to_string()),
            b'I'=>three_prot.push("Ile".to_string()),
            b'M'=>three_prot.push("Met".to_string()),
            b'T'=>three_prot.push("Thr".to_string()),
            b'N'=>three_prot.push("Asn".to_string()),
            b'K'=>three_prot.push("Lys".to_string()),
            b'V'=>three_prot.push("Val".to_string()),
            b'A'=>three_prot.push("Ala".to_string()),
            b'D'=>three_prot.push("Asp".to_string()),
            b'E'=>three_prot.push("Glu".to_string()),
            b'G'=>three_prot.push("Gly".to_string()),
            _=>three_prot.push("Wrong".to_string()),
        }
    }
    three_prot
}


///Converts a normal string, like FASTA format without the first line, to the Text vector format
pub fn raw_strand_to_textformat(strand : String) -> Text {
    let mut vec_strand : Text = vec![];
    for i in strand.into_bytes() {
        vec_strand.push(i);
    }
    vec_strand
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::ops::Deref;

    fn print_sequence<Item: Deref<Target = u8>, T: IntoIterator<Item = Item>>(sequence: T) {
        for c in sequence {
            println!("{}", *c);
        }
    }

    #[test]
    fn printing_strand() {
        let seq_ex = b"ACGT";
        let seq_vec = &vec![b'A',b'C'];
        print_sequence(seq_ex.iter().step_by(1));
        print_sequence(seq_vec);
    }

    #[test]
    fn transcript_simple_dna(){
        let seq_vec : Text= vec![b'A',b'C',b'T'];
        let new_arn = transcription(seq_vec);
        assert_eq!(new_arn,vec![b'A',b'C',b'U']);
        print_sequence(&new_arn);
    }

    #[test]
    fn traduction_one_codon_protein(){
        let seq_arn : Text = vec![b'U',b'U',b'U'];
        let new_prot : Text = traduction(seq_arn);
        assert_eq!(new_prot,vec![b'F']);
        print_sequence(&new_prot);
    }

    #[test]
    fn uni_to_three_protein(){
        let prot_ex : Text = vec![b'F',b'T'];
        let prot_result = protein_to_three(prot_ex);
        assert_eq!(prot_result,vec!["Phe".to_string(),"Thr".to_string()]);
    }

    #[test]
    fn string_to_textformat(){
        let str_initial = String::from("ATGCCGGTTAGCAACGT");
        let str_result = raw_strand_to_textformat(str_initial);
        let str_expect : Text = vec![b'A',b'T',b'G',b'C',b'C',b'G',b'G',b'T',b'T',b'A',b'G',b'C',b'A',b'A',b'C',b'G',b'T'];
        assert_eq!(str_expect,str_result);
    }
}

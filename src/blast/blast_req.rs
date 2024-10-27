//! Request with http method to the NCBI API 
//!  
//!  ```
//!  use bio::blast::blast_req;
//!  let dna_string : String = String::from("ATGGTAAGTTATTCTTCTAGTTGATATATTCTTACTCTTTCTTTCTTACGTAACAACTGATCAAGTGTTAATATATTGTATTTAATCAGCAAGGCCAGTGGATAGCCGCAAGAGACCTTTCTATTACATGGGTGGACAATCCTCAGTACTGGACATGGAAAACTGTTGATCCTAAGTCTGTCTCTCTCTAACTATCTCTCTCTAAATTTCATTAAAAAGAGAAATATATTTTAATCCGATCTTAACATTTTGTACTACAGTATTGAAGTGGCAGAGCTTCTTAGGGTAGCTTGGCTTGACATTTATGGAAAGATCGAGACAAAAAATCTTATTCGAAAGACTAGTTATGCTGTATATTTAGTGTTCAAGTTAACAGATAACCCTCGTGAACTTGAACGAGCCACAGCGTCGCTAAGATTTGTGAACGAAGTGGCGGAGGGCGCTGGCATTGAGGGTACCACTGTTTTCATCTCGAAGAAAAAGGAATTACCAGGAGAACTTGGCCGGTTCCCACATCTCCGAAGTGATGGCTGGTTAGAAACCAAGCTTGGTGAGTTTTTCAACAACTTAGGAGAGGATGGTGAAGTCGAAATGAGGTTGATGGAAATCAATGACAAAACTTGGAAATCTGGCATCATTGTTAAGGGCTTCGACATTCGTCCAAACTAA");
//!  let dna_vector : Text = blast::raw_strand_to_textformat(dna_string);
//! 
//!  let fasta_format : String = blast_req::create_fasta_from_dna(dna_vector); //Can be converted 'manually'
//! 
//!  let blastn_results : String = blast_req::blastn_req(dna_vector);
//! 
//!  println!("{:#?}",blastn_results);
//!   
//!  ```
//!  


use std::error::Error;
use crate::utils::Text;
use reqwest::{self};
use serde_json::Value;
use std::fmt;

#[derive(Debug)]
struct MyError(String);

///Error implemented for http request methology required for the blastn request
impl fmt::Display for MyError{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "There is an error: \n {}", self.0)
    }
}

impl Error for MyError {}

static URL_NCBI: &'static str = "https://blast.ncbi.nlm.nih.gov/Blast.cgi";

///Transform the nucleotide strand in Text format to the FASTA format on string supported by NCBI site
pub fn create_fasta_from_dna(dna_str : Text) -> String {
    let mut dna_fasta = String::from(">Search dna from rust-bio \n");
    let mut count = 0;
    for i in dna_str{
        let nc = i as char;
        if count < 69{
            dna_fasta.push(nc);
            count=count+1;
        }
        else {
            count = 0;
            dna_fasta.push('\n');
            dna_fasta.push(nc);
        }
    }

    dna_fasta.to_owned()
}

///Load a search on the blastn platform and after a positive response is given, the text format of the main result is returned
pub async fn send_http_req(dna_fasta : String) -> Result<String,Box<dyn Error>> {
    if !dna_fasta.contains("\n") || !dna_fasta.contains('>') {
        return Err(Box::new(MyError("The string must be on FASTA format".into())));
    }
    let client_ncbi = reqwest::Client::new();
    let params = [
        (String::from("QUERY"),dna_fasta),
        (String::from("DATABASE"),String::from("nt")),
        (String::from("PROGRAM"),String::from("blastn")),
        (String::from("CMD"),String::from("Put"))
        ];
    //let mut url = reqwest::Url::parse_with_params(URL_NCBI, &params)?;
    let resp = client_ncbi.get(URL_NCBI).query(&params)
        .send().await.or(Err(Box::new(MyError("No response was recieved from client".into()))))?;
    let resp_text = resp.text().await
        .or(Err(Box::new(MyError("Response was not able to convert to text".into()))))?;
    //println!("{:#?}", resp_text);

    let index_rid = resp_text.find("RID = ").unwrap();
    let start_rid = index_rid+6;
    let iend_rid= resp_text.find("RTOE = ").unwrap();
    let end_rid = iend_rid-5;
    let rid = &resp_text[start_rid..end_rid];

    if rid.len() != 11{
        return Err(Box::new(MyError("Length of rid is wrong {rid}".into())));
    }
    let params_result = [
        ("CMD","Get"),
        ("RID",rid)
    ];
    let mut status = String::from("WAITI");
    let mut search_resp : reqwest::Response;

    while status != "READY"{
        search_resp = client_ncbi.get(URL_NCBI)
            .query(&params_result).send().await
            .or(Err(Box::new(MyError("No response was recieved from client with the ID while waiting".into()))))?;
        let resp_search_text = search_resp.text().await
            .or(Err(Box::new(MyError("Response was not able to convert to text while waiting".into()))))?;
        let index_status = resp_search_text.find("Status=").unwrap() +7;
        let max_index = index_status + 5;
        status.clear();
        status = (&resp_search_text[index_status..max_index]).to_string();
        if status == "UNKWO"{
            return Err(Box::new(MyError("The ID saved was not found on the server".into())));
        }
        else{continue;}
        
    }

    let json_params_result = [
        ("CMD","Get"),
        ("RID",rid),
        ("FORMAT_TYPE","Text"),
        ("ALIGNMENTS","1"),
        ("DESCRIPTIONS","1"),
    ];
    search_resp = client_ncbi.get(URL_NCBI)
        .query(&json_params_result).send().await
        .or(Err(Box::new(MyError("Error on the response requested as JSON format".into()))))?;
    let final_text = search_resp.text().await
        .or(Err(Box::new(MyError("Error converting response JSON response to text".into()))))?;

    /*let parsed = read_json(&json_text);
    let list_res = parsed["hits"][0].as_str().clone();
    let blast_result = list_res.unwrap().to_string();*/
    println!("{:#?}", final_text);
    let start_after_box = final_text.find("\nBLASTN").unwrap() +1;
    let search_result_clear = final_text[start_after_box..].to_string();


    if search_result_clear.is_empty() {
        Err(Box::new(MyError("No match found".into())))?
    }else{
    Ok(search_result_clear)}
}

///Get a json alike format from a json text
fn _read_json(raw_json : &str) -> Value {
    let parsed: Value = serde_json::from_str(raw_json).unwrap();
    parsed
}

///Accepts a Text vector of nucleotides which is used to request a search on NCBI blastNucleotide. This returns a string with the text format of the main result given by the most alike nucleotide strand on the database.
pub async fn blastn_req(dna_str : Text) -> String {
    let dna_fasta = create_fasta_from_dna(dna_str);
    let results_blastn = send_http_req(dna_fasta).await;
    match results_blastn {
        Ok(n) => {return n},
        Err(e) => {
            panic!("Error found {}", e)
        },
    }
    
}

#[cfg(test)]
#[allow(non_snake_case)]
mod tests {
    use crate::blast::blast_req::send_http_req;
    use crate::blast::Text;

    use super::{create_fasta_from_dna, blastn_req};

    #[test]
    fn print_html_response() {
        let fasta_example = String::from(">Example
        ATGGTAAGTTATTCTTCTAGTTGATATATTCTTACTCTTTCTTTCTTACGTAACAACTGATCAAGTGTTA
        ATATATTGTATTTAATCAGCAAGGCCAGTGGATAGCCGCAAGAGACCTTTCTATTACATGGGTGGACAAT
        CCTCAGTACTGGACATGGAAAACTGTTGATCCTAAGTCTGTCTCTCTCTAACTATCTCTCTCTAAATTTC
        ATTAAAAAGAGAAATATATTTTAATCCGATCTTAACATTTTGTACTACAGTATTGAAGTGGCAGAGCTTC
        TTAGGGTAGCTTGGCTTGACATTTATGGAAAGATCGAGACAAAAAATCTTATTCGAAAGACTAGTTATGC
        TGTATATTTAGTGTTCAAGTTAACAGATAACCCTCGTGAACTTGAACGAGCCACAGCGTCGCTAAGATTT
        GTGAACGAAGTGGCGGAGGGCGCTGGCATTGAGGGTACCACTGTTTTCATCTCGAAGAAAAAGGAATTAC
        CAGGAGAACTTGGCCGGTTCCCACATCTCCGAAGTGATGGCTGGTTAGAAACCAAGCTTGGTGAGTTTTT
        CAACAACTTAGGAGAGGATGGTGAAGTCGAAATGAGGTTGATGGAAATCAATGACAAAACTTGGAAATCT
        GGCATCATTGTTAAGGGCTTCGACATTCGTCCAAACTAA");
        
        let _result = send_http_req(fasta_example);
    }

    #[test]
    fn dna_to_fasta(){
        let dna_example :Text = vec![b'A',b'T',b'G',b'G',b'T',b'A',b'A',b'G',b'T',b'T',b'A',b'T',b'T',b'C',b'T',b'T',b'C',b'T',b'A',b'G',b'T',b'T',b'G',b'A',b'T',b'A',b'T',b'A',b'T',b'T',b'C',b'T',b'T',b'A',b'C',b'T',b'C',b'T',b'T',b'T',b'C',b'T',b'T',b'T',b'C',b'T',b'T',b'A',b'C',b'G',b'T',b'A',b'A',b'C',b'A',b'A',b'C',b'T',b'G',b'A',b'T',b'C',b'A',b'A',b'G',b'T',b'G',b'T',b'T',b'A',
        b'A',b'T',b'A',b'T',b'A',b'T',b'T',b'G',b'T',b'A',b'T',b'T',b'T',b'A',b'A',b'T',b'C',b'A',b'G',b'C',b'A',b'A',b'G',b'G',b'C',b'C',b'A',b'G',b'T',b'G',b'G',b'A',b'T',b'A',b'G',b'C',b'C',b'G',b'C',b'A',b'A',b'G',b'A',b'G',b'A',b'C',b'C',b'T',b'T',b'T',b'C',b'T',b'A',b'T',b'T',b'A',b'C',b'A',b'T',b'G',b'G',b'G',b'T',b'G',b'G',b'A',b'C',b'A',b'A',b'T',
        b'C',b'C',b'T',b'C',b'A',b'G',b'T',b'A',b'C',b'T',b'G',b'G',b'A',b'C',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'A',b'C',b'T',b'G',b'T',b'T',b'G',b'A',b'T',b'C',b'C',b'T',b'A',b'A',b'G',b'T',b'C',b'T',b'G',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'A',b'A',b'C',b'T',b'A',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'A',b'A',b'A',b'T',b'T',b'T',b'C',
        b'A',b'T',b'T',b'A',b'A',b'A',b'A',b'A',b'G',b'A',b'G',b'A',b'A',b'A',b'T',b'A',b'T',b'A',b'T',b'T',b'T',b'T',b'A',b'A',b'T',b'C',b'C',b'G',b'A',b'T',b'C',b'T',b'T',b'A',b'A',b'C',b'A',b'T',b'T',b'T',b'T',b'G',b'T',b'A',b'C',b'T',b'A',b'C',b'A',b'G',b'T',b'A',b'T',b'T',b'G',b'A',b'A',b'G',b'T',b'G',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'T',b'C',
        b'T',b'T',b'A',b'G',b'G',b'G',b'T',b'A',b'G',b'C',b'T',b'T',b'G',b'G',b'C',b'T',b'T',b'G',b'A',b'C',b'A',b'T',b'T',b'T',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'G',b'A',b'T',b'C',b'G',b'A',b'G',b'A',b'C',b'A',b'A',b'A',b'A',b'A',b'A',b'T',b'C',b'T',b'T',b'A',b'T',b'T',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'C',b'T',b'A',b'G',b'T',b'T',b'A',b'T',b'G',b'C',
        b'T',b'G',b'T',b'A',b'T',b'A',b'T',b'T',b'T',b'A',b'G',b'T',b'G',b'T',b'T',b'C',b'A',b'A',b'G',b'T',b'T',b'A',b'A',b'C',b'A',b'G',b'A',b'T',b'A',b'A',b'C',b'C',b'C',b'T',b'C',b'G',b'T',b'G',b'A',b'A',b'C',b'T',b'T',b'G',b'A',b'A',b'C',b'G',b'A',b'G',b'C',b'C',b'A',b'C',b'A',b'G',b'C',b'G',b'T',b'C',b'G',b'C',b'T',b'A',b'A',b'G',b'A',b'T',b'T',b'T',
        b'G',b'T',b'G',b'A',b'A',b'C',b'G',b'A',b'A',b'G',b'T',b'G',b'G',b'C',b'G',b'G',b'A',b'G',b'G',b'G',b'C',b'G',b'C',b'T',b'G',b'G',b'C',b'A',b'T',b'T',b'G',b'A',b'G',b'G',b'G',b'T',b'A',b'C',b'C',b'A',b'C',b'T',b'G',b'T',b'T',b'T',b'T',b'C',b'A',b'T',b'C',b'T',b'C',b'G',b'A',b'A',b'G',b'A',b'A',b'A',b'A',b'A',b'G',b'G',b'A',b'A',b'T',b'T',b'A',b'C',
        b'C',b'A',b'G',b'G',b'A',b'G',b'A',b'A',b'C',b'T',b'T',b'G',b'G',b'C',b'C',b'G',b'G',b'T',b'T',b'C',b'C',b'C',b'A',b'C',b'A',b'T',b'C',b'T',b'C',b'C',b'G',b'A',b'A',b'G',b'T',b'G',b'A',b'T',b'G',b'G',b'C',b'T',b'G',b'G',b'T',b'T',b'A',b'G',b'A',b'A',b'A',b'C',b'C',b'A',b'A',b'G',b'C',b'T',b'T',b'G',b'G',b'T',b'G',b'A',b'G',b'T',b'T',b'T',b'T',b'T',
        b'C',b'A',b'A',b'C',b'A',b'A',b'C',b'T',b'T',b'A',b'G',b'G',b'A',b'G',b'A',b'G',b'G',b'A',b'T',b'G',b'G',b'T',b'G',b'A',b'A',b'G',b'T',b'C',b'G',b'A',b'A',b'A',b'T',b'G',b'A',b'G',b'G',b'T',b'T',b'G',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'T',b'C',b'A',b'A',b'T',b'G',b'A',b'C',b'A',b'A',b'A',b'A',b'C',b'T',b'T',b'G',b'G',b'A',b'A',b'A',b'T',b'C',b'T',
        b'G',b'G',b'C',b'A',b'T',b'C',b'A',b'T',b'T',b'G',b'T',b'T',b'A',b'A',b'G',b'G',b'G',b'C',b'T',b'T',b'C',b'G',b'A',b'C',b'A',b'T',b'T',b'C',b'G',b'T',b'C',b'C',b'A',b'A',b'A',b'C',b'T',b'A',b'A'];
        
        let result_fasta = create_fasta_from_dna(dna_example);
        let fasta_expected = String::from(">Search dna from rust-bio \nATGGTAAGTTATTCTTCTAGTTGATATATTCTTACTCTTTCTTTCTTACGTAACAACTGATCAAGTGTT
AATATATTGTATTTAATCAGCAAGGCCAGTGGATAGCCGCAAGAGACCTTTCTATTACATGGGTGGACAA
TCCTCAGTACTGGACATGGAAAACTGTTGATCCTAAGTCTGTCTCTCTCTAACTATCTCTCTCTAAATTT
CATTAAAAAGAGAAATATATTTTAATCCGATCTTAACATTTTGTACTACAGTATTGAAGTGGCAGAGCTT
CTTAGGGTAGCTTGGCTTGACATTTATGGAAAGATCGAGACAAAAAATCTTATTCGAAAGACTAGTTATG
CTGTATATTTAGTGTTCAAGTTAACAGATAACCCTCGTGAACTTGAACGAGCCACAGCGTCGCTAAGATT
TGTGAACGAAGTGGCGGAGGGCGCTGGCATTGAGGGTACCACTGTTTTCATCTCGAAGAAAAAGGAATTA
CCAGGAGAACTTGGCCGGTTCCCACATCTCCGAAGTGATGGCTGGTTAGAAACCAAGCTTGGTGAGTTTT
TCAACAACTTAGGAGAGGATGGTGAAGTCGAAATGAGGTTGATGGAAATCAATGACAAAACTTGGAAATC
TGGCATCATTGTTAAGGGCTTCGACATTCGTCCAAACTAA");
        assert_eq!(result_fasta,fasta_expected);
    }

    #[actix_rt::test]
    async fn send_basic_fasta() {
        let dna_text :Text = vec![b'A',b'T',b'G',b'G',b'T',b'A',b'A',b'G',b'T',b'T',b'A',b'T',b'T',b'C',b'T',b'T',b'C',b'T',b'A',b'G',b'T',b'T',b'G',b'A',b'T',b'A',b'T',b'A',b'T',b'T',b'C',b'T',b'T',b'A',b'C',b'T',b'C',b'T',b'T',b'T',b'C',b'T',b'T',b'T',b'C',b'T',b'T',b'A',b'C',b'G',b'T',b'A',b'A',b'C',b'A',b'A',b'C',b'T',b'G',b'A',b'T',b'C',b'A',b'A',b'G',b'T',b'G',b'T',b'T',b'A',
        b'A',b'T',b'A',b'T',b'A',b'T',b'T',b'G',b'T',b'A',b'T',b'T',b'T',b'A',b'A',b'T',b'C',b'A',b'G',b'C',b'A',b'A',b'G',b'G',b'C',b'C',b'A',b'G',b'T',b'G',b'G',b'A',b'T',b'A',b'G',b'C',b'C',b'G',b'C',b'A',b'A',b'G',b'A',b'G',b'A',b'C',b'C',b'T',b'T',b'T',b'C',b'T',b'A',b'T',b'T',b'A',b'C',b'A',b'T',b'G',b'G',b'G',b'T',b'G',b'G',b'A',b'C',b'A',b'A',b'T',
        b'C',b'C',b'T',b'C',b'A',b'G',b'T',b'A',b'C',b'T',b'G',b'G',b'A',b'C',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'A',b'C',b'T',b'G',b'T',b'T',b'G',b'A',b'T',b'C',b'C',b'T',b'A',b'A',b'G',b'T',b'C',b'T',b'G',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'A',b'A',b'C',b'T',b'A',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'C',b'T',b'A',b'A',b'A',b'T',b'T',b'T',b'C',
        b'A',b'T',b'T',b'A',b'A',b'A',b'A',b'A',b'G',b'A',b'G',b'A',b'A',b'A',b'T',b'A',b'T',b'A',b'T',b'T',b'T',b'T',b'A',b'A',b'T',b'C',b'C',b'G',b'A',b'T',b'C',b'T',b'T',b'A',b'A',b'C',b'A',b'T',b'T',b'T',b'T',b'G',b'T',b'A',b'C',b'T',b'A',b'C',b'A',b'G',b'T',b'A',b'T',b'T',b'G',b'A',b'A',b'G',b'T',b'G',b'G',b'C',b'A',b'G',b'A',b'G',b'C',b'T',b'T',b'C',
        b'T',b'T',b'A',b'G',b'G',b'G',b'T',b'A',b'G',b'C',b'T',b'T',b'G',b'G',b'C',b'T',b'T',b'G',b'A',b'C',b'A',b'T',b'T',b'T',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'G',b'A',b'T',b'C',b'G',b'A',b'G',b'A',b'C',b'A',b'A',b'A',b'A',b'A',b'A',b'T',b'C',b'T',b'T',b'A',b'T',b'T',b'C',b'G',b'A',b'A',b'A',b'G',b'A',b'C',b'T',b'A',b'G',b'T',b'T',b'A',b'T',b'G',b'C',
        b'T',b'G',b'T',b'A',b'T',b'A',b'T',b'T',b'T',b'A',b'G',b'T',b'G',b'T',b'T',b'C',b'A',b'A',b'G',b'T',b'T',b'A',b'A',b'C',b'A',b'G',b'A',b'T',b'A',b'A',b'C',b'C',b'C',b'T',b'C',b'G',b'T',b'G',b'A',b'A',b'C',b'T',b'T',b'G',b'A',b'A',b'C',b'G',b'A',b'G',b'C',b'C',b'A',b'C',b'A',b'G',b'C',b'G',b'T',b'C',b'G',b'C',b'T',b'A',b'A',b'G',b'A',b'T',b'T',b'T',
        b'G',b'T',b'G',b'A',b'A',b'C',b'G',b'A',b'A',b'G',b'T',b'G',b'G',b'C',b'G',b'G',b'A',b'G',b'G',b'G',b'C',b'G',b'C',b'T',b'G',b'G',b'C',b'A',b'T',b'T',b'G',b'A',b'G',b'G',b'G',b'T',b'A',b'C',b'C',b'A',b'C',b'T',b'G',b'T',b'T',b'T',b'T',b'C',b'A',b'T',b'C',b'T',b'C',b'G',b'A',b'A',b'G',b'A',b'A',b'A',b'A',b'A',b'G',b'G',b'A',b'A',b'T',b'T',b'A',b'C',
        b'C',b'A',b'G',b'G',b'A',b'G',b'A',b'A',b'C',b'T',b'T',b'G',b'G',b'C',b'C',b'G',b'G',b'T',b'T',b'C',b'C',b'C',b'A',b'C',b'A',b'T',b'C',b'T',b'C',b'C',b'G',b'A',b'A',b'G',b'T',b'G',b'A',b'T',b'G',b'G',b'C',b'T',b'G',b'G',b'T',b'T',b'A',b'G',b'A',b'A',b'A',b'C',b'C',b'A',b'A',b'G',b'C',b'T',b'T',b'G',b'G',b'T',b'G',b'A',b'G',b'T',b'T',b'T',b'T',b'T',
        b'C',b'A',b'A',b'C',b'A',b'A',b'C',b'T',b'T',b'A',b'G',b'G',b'A',b'G',b'A',b'G',b'G',b'A',b'T',b'G',b'G',b'T',b'G',b'A',b'A',b'G',b'T',b'C',b'G',b'A',b'A',b'A',b'T',b'G',b'A',b'G',b'G',b'T',b'T',b'G',b'A',b'T',b'G',b'G',b'A',b'A',b'A',b'T',b'C',b'A',b'A',b'T',b'G',b'A',b'C',b'A',b'A',b'A',b'A',b'C',b'T',b'T',b'G',b'G',b'A',b'A',b'A',b'T',b'C',b'T',
        b'G',b'G',b'C',b'A',b'T',b'C',b'A',b'T',b'T',b'G',b'T',b'T',b'A',b'A',b'G',b'G',b'G',b'C',b'T',b'T',b'C',b'G',b'A',b'C',b'A',b'T',b'T',b'C',b'G',b'T',b'C',b'C',b'A',b'A',b'A',b'C',b'T',b'A',b'A'];
        
        let result = blastn_req(dna_text).await;
        let result_expected = String::from("BLASTN 2.14.1+\nReference: Stephen F. Altschul, Thomas L. Madden, Alejandro\nA. Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and\nDavid J. Lipman (1997), \"Gapped BLAST and PSI-BLAST: a new\ngeneration of protein database search programs\", Nucleic\nAcids Res. 25:3389-3402.\n\n\nRID: ABM23074013\n\n\nDatabase: Nucleotide collection (nt)\n           95,707,385 sequences; 1,224,912,062,916 total letters\nQuery= Search dna from rust-bio\n\nLength=669\n\n\n                                                                   Score     E     Max\nSequences producing significant alignments:                       (Bits)  Value  Ident\n\nDQ155404.1 Nicotiana tomentosa Tom gene, complete cds              1207    0.0    100%        \n\nALIGNMENTS\n>DQ155404.1 Nicotiana tomentosa Tom gene, complete cds\nLength=669\n\n Score = 1207 bits (1338),  Expect = 0.0\n Identities = 669/669 (100%), Gaps = 0/669 (0%)\n Strand=Plus/Plus\n\nQuery  1    ATGGTAAGTTATTCTTCTAGTTGATATATTCTTACTCTTTCTTTCTTACGTAACAACTGA  60\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  1    ATGGTAAGTTATTCTTCTAGTTGATATATTCTTACTCTTTCTTTCTTACGTAACAACTGA  60\n\nQuery  61   TCAAGTGTTAATATATTGTATTTAATCAGCAAGGCCAGTGGATAGCCGCAAGAGACCTTT  120\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  61   TCAAGTGTTAATATATTGTATTTAATCAGCAAGGCCAGTGGATAGCCGCAAGAGACCTTT  120\n\nQuery  121  CTATTACATGGGTGGACAATCCTCAGTACTGGACATGGAAAACTGTTGATCCTAAGTCTG  180\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  121  CTATTACATGGGTGGACAATCCTCAGTACTGGACATGGAAAACTGTTGATCCTAAGTCTG  180\n\nQuery  181  tctctctctaactatctctctctAAATTTCATTAAAAAGAGAAATATATTTTAATCCGAT  240\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  181  TCTCTCTCTAACTATCTCTCTCTAAATTTCATTAAAAAGAGAAATATATTTTAATCCGAT  240\n\nQuery  241  CTTAACATTTTGTACTACAGTATTGAAGTGGCAGAGCTTCTTAGGGTAGCTTGGCTTGAC  300\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  241  CTTAACATTTTGTACTACAGTATTGAAGTGGCAGAGCTTCTTAGGGTAGCTTGGCTTGAC  300\n\nQuery  301  ATTTATGGAAAGATCGAGACAAAAAATCTTATTCGAAAGACTAGTTATGCTGTATATTTA  360\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  301  ATTTATGGAAAGATCGAGACAAAAAATCTTATTCGAAAGACTAGTTATGCTGTATATTTA  360\n\nQuery  361  GTGTTCAAGTTAACAGATAACCCTCGTGAACTTGAACGAGCCACAGCGTCGCTAAGATTT  420\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  361  GTGTTCAAGTTAACAGATAACCCTCGTGAACTTGAACGAGCCACAGCGTCGCTAAGATTT  420\n\nQuery  421  GTGAACGAAGTGGCGGAGGGCGCTGGCATTGAGGGTACCACTGTTTTCATCTCGAAGAAA  480\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  421  GTGAACGAAGTGGCGGAGGGCGCTGGCATTGAGGGTACCACTGTTTTCATCTCGAAGAAA  480\n\nQuery  481  AAGGAATTACCAGGAGAACTTGGCCGGTTCCCACATCTCCGAAGTGATGGCTGGTTAGAA  540\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  481  AAGGAATTACCAGGAGAACTTGGCCGGTTCCCACATCTCCGAAGTGATGGCTGGTTAGAA  540\n\nQuery  541  ACCAAGCTTGGTGAGTTTTTCAACAACTTAGGAGAGGATGGTGAAGTCGAAATGAGGTTG  600\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  541  ACCAAGCTTGGTGAGTTTTTCAACAACTTAGGAGAGGATGGTGAAGTCGAAATGAGGTTG  600\n\nQuery  601  ATGGAAATCAATGACAAAACTTGGAAATCTGGCATCATTGTTAAGGGCTTCGACATTCGT  660\n            ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||\nSbjct  601  ATGGAAATCAATGACAAAACTTGGAAATCTGGCATCATTGTTAAGGGCTTCGACATTCGT  660\n\nQuery  661  CCAAACTAA  669\n            |||||||||\nSbjct  661  CCAAACTAA  669\n\n\n  Database: Nucleotide collection (nt)\n    Posted date:  Jul 1, 2023  4:29 PM\n  Number of letters in database: 1,224,912,062,916\n  Number of sequences in database:  95,707,385\n\nLambda      K        H\n   0.634    0.408    0.912 \nGapped\nLambda      K        H\n   0.625    0.410    0.780 \nMatrix: blastn matrix:2 -3\nGap Penalties: Existence: 5, Extension: 2\nNumber of Sequences: 95707385\nNumber of Hits to DB: 45866596\nNumber of extensions: 101342\nNumber of successful extensions: 101342\nNumber of sequences better than 10: 21\nNumber of HSP's better than 10 without gapping: 0\nNumber of HSP's gapped: 101328\nNumber of HSP's successfully gapped: 29\nLength of query: 669\nLength of database: 1224912062916\nLength adjustment: 40\nEffective length of query: 629\nEffective length of database: 1221083767516\nEffective search space: 768061689767564\nEffective search space used: 768061689767564\nA: 0\nX1: 22 (20.1 bits)\nX2: 33 (29.8 bits)\nX3: 110 (99.2 bits)\nS1: 28 (26.5 bits)\nS2: 50 (46.4 bits)\n\n\n");
        let res_first_20 = &result[..20];
        let exp_first_20 = &result_expected[..20];
        assert_eq!(res_first_20,exp_first_20);

        let size_result = result.len();
        assert!(size_result >= 200);

    }
}
/*
1 Func e volver adn a fasta DONE
2 generar el request
3 blast req principal publico
4 aplicarlo a blast p tamb
 */
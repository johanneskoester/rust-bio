// Copyright 2019 Magnus Manske
// Licensed under the MIT license (http://opensource.org/licenses/MIT)
// This file may not be copied, modified, or distributed
// except according to those terms.

//! GAF2 format reading and writing.
//!
//! GAF2.1 definition : https://geneontology.github.io/docs/go-annotation-file-gaf-format-21/
//!
//! # Example
//!
//! ```
//! use std::io;
//! use bio::io::gaf;
//! let reader = gaf::Reader::new(io::stdin(), gaf::GafType::GAF2);
//! ```

use itertools::Itertools;
use multimap::MultiMap;
use regex::Regex;
use std::convert::AsRef;
use std::fmt;
use std::fs;
use std::io;
use std::path::Path;

use csv;

/// `GafType`
///
/// We have one format in the GAF family.
/// The change is in the last field of GAF.
/// For each type we have key value separator and field separator
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum GafType {
    GAF1,
    GAF2,
}

impl GafType {
    #[inline]
    fn separator(&self) -> (u8, u8, u8) {
        match *self {
            GafType::GAF1 => (b':', b'|', 0u8),
            GafType::GAF2 => (b':', b'|', 0u8),
        }
    }

    // Used in writer. Should be GafType::GAF1 instead of string in comparison there, but can't quite make it work
    #[inline]
    pub fn name(&self) -> String {
        match *self {
            GafType::GAF1 => "GAF1".to_string(),
            GafType::GAF2 => "GAF2".to_string(),
        }
    }
}

/// A GAF reader.
#[derive(Debug)]
pub struct Reader<R: io::Read> {
    inner: csv::Reader<R>,
    gaf_type: GafType,
}

impl Reader<fs::File> {
    /// Read GAF from given file path in given format.
    pub fn from_file<P: AsRef<Path>>(path: P, fileformat: GafType) -> io::Result<Self> {
        fs::File::open(path).map(|f| Reader::new(f, fileformat))
    }
}

impl<R: io::Read> Reader<R> {
    /// Create a new GAF reader given an instance of `io::Read`, in given format.
    pub fn new(reader: R, fileformat: GafType) -> Self {
        Reader {
            inner: csv::ReaderBuilder::new()
                .delimiter(b'\t')
                .has_headers(false)
                .comment(Some(b'!'))
                .from_reader(reader),
            gaf_type: fileformat,
        }
    }

    /// Iterate over all records.
    pub fn records(&mut self) -> Records<R> {
        let (delim, term, vdelim) = self.gaf_type.separator();
        let r1 = format!(
            r" *(?P<key>[^{delim}\{term}\t]+){delim}(?P<value>[^{delim}\{term}\t]+)\{term}?",
            //r" *(?P<key>[^:\|\t]+):(?P<value>[^:\|\t]+)\|?",
            delim = delim as char,
            term = term as char
        );
        let db_ref_re = Regex::new(&r1).unwrap();

        let r2 = format!(r" *(?P<value>[^{term}\t]+)\{term}?", term = term as char);
        let default_re = Regex::new(&r2).unwrap();

        Records {
            inner: self.inner.deserialize(),
            db_ref_re,
            default_re,
            value_delim: vdelim as char,
        }
    }
}

/// A GAF record.
pub struct Records<'a, R: 'a + io::Read> {
    inner: csv::DeserializeRecordsIter<'a, R, Vec<String>>,
    db_ref_re: Regex,
    default_re: Regex,
    value_delim: char,
}

impl<'a, R: io::Read> Iterator for Records<'a, R> {
    type Item = csv::Result<Record>;

    fn next(&mut self) -> Option<csv::Result<Record>> {
        self.inner
            .next()
            .map(|res| res.map(|v| self.from_csv_row(v)))
    }
}

impl<'a, R: io::Read> Records<'a, R> {
    fn from_csv_row(&mut self, v: Vec<String>) -> Record {
        let trim_quotes = |s: &str| s.trim_matches('\'').trim_matches('"').to_owned();

        let db = v[0].clone();
        let db_object_id = v[1].clone();
        let db_object_symbol = v[2].clone();
        let raw_qualifier = v[3].clone();
        let go_id = v[4].clone();
        let raw_db_ref = v[5].clone();
        let evidence_code = v[6].clone();
        let raw_with_from = v[7].clone();
        let aspect = v[8].clone();
        let db_object_name = v[9].clone();
        let raw_db_object_synonym = v[10].clone();
        let db_object_type = v[11].clone();
        let raw_taxon = v[12].clone();
        let date = v[13].clone();
        let assigned_by = v[14].clone();
        let raw_annotation_extension = if v.len() >= 16 {
            v[15].clone()
        } else {
            String::new()
        };
        let gene_product_form_id: Option<String> = if v.len() >= 17 {
            Some(v[16].clone())
        } else {
            None
        };

        let mut db_ref = MultiMap::new();
        for caps in self.db_ref_re.captures_iter(&raw_db_ref.to_string()) {
            for value in caps["value"].split(self.value_delim) {
                db_ref.insert(trim_quotes(&caps["key"]), trim_quotes(value));
            }
        }

        let mut qualifier = Vec::<String>::new();
        for caps in self.default_re.captures_iter(&raw_qualifier.to_string()) {
            qualifier.push(trim_quotes(&caps["value"]));
        }

        let mut with_from = Vec::<String>::new();
        for caps in self.default_re.captures_iter(&raw_with_from.to_string()) {
            with_from.push(trim_quotes(&caps["value"]));
        }

        let mut db_object_synonym = Vec::<String>::new();
        for caps in self
            .default_re
            .captures_iter(&raw_db_object_synonym.to_string())
        {
            db_object_synonym.push(trim_quotes(&caps["value"]));
        }

        let mut taxon = Vec::<String>::new();
        for caps in self.default_re.captures_iter(&raw_taxon.to_string()) {
            taxon.push(trim_quotes(&caps["value"]));
        }

        let mut annotation_extension: Option<Vec<String>> = None;
        for caps in self
            .default_re
            .captures_iter(&raw_annotation_extension.to_string())
        {
            if annotation_extension.is_some() {
                match annotation_extension.as_mut() {
                    Some(ae) => {
                        ae.push(trim_quotes(&caps["value"]));
                    }
                    None => {}
                }
            } else {
                annotation_extension = Some(vec![trim_quotes(&caps["value"])]);
            }
        }

        Record {
            db,
            db_object_id,
            db_object_symbol,
            qualifier,
            go_id,
            db_ref,
            evidence_code,
            with_from,
            aspect,
            db_object_name,
            db_object_synonym,
            db_object_type,
            taxon,
            date,
            assigned_by,
            annotation_extension,
            gene_product_form_id,
        }
    }
}

/// A GAF writer.
#[derive(Debug)]
pub struct Writer<W: io::Write> {
    inner: csv::Writer<W>,
    delimiter: char,
    terminator: String,
    fileformat: GafType,
}

impl Writer<fs::File> {
    /// Write to a given file path in given format.
    pub fn to_file<P: AsRef<Path>>(path: P, fileformat: GafType) -> io::Result<Self> {
        fs::File::create(path).map(|f| Writer::new(f, fileformat))
    }
}

impl<W: io::Write> Writer<W> {
    /// Write to a given writer.
    pub fn new(writer: W, fileformat: GafType) -> Self {
        let (delim, termi, _) = fileformat.separator();

        Writer {
            inner: csv::WriterBuilder::new()
                .delimiter(b'\t')
                .flexible(true)
                .from_writer(writer),
            delimiter: delim as char,
            terminator: String::from_utf8(vec![termi]).unwrap(),
            fileformat: fileformat,
        }
    }

    /// Write a given GAF record.
    pub fn write(&mut self, record: &Record) -> csv::Result<()> {
        let v = record.into_csv(
            &self.fileformat,
            self.delimiter.to_string(),
            self.terminator.to_string(),
        );
        self.inner.serialize(v)
    }
}

/// A GAF record
#[derive(Debug, Default, Serialize, Deserialize, Clone)]
pub struct Record {
    db: String,
    db_object_id: String,
    db_object_symbol: String,
    qualifier: Vec<String>,
    go_id: String,
    db_ref: MultiMap<String, String>,
    evidence_code: String,
    with_from: Vec<String>,
    aspect: String,
    db_object_name: String,
    db_object_synonym: Vec<String>,
    db_object_type: String,
    taxon: Vec<String>,
    date: String,
    assigned_by: String,
    annotation_extension: Option<Vec<String>>,
    gene_product_form_id: Option<String>,
}

/// Implementing fmt::Display for easier output. Assumes GAF2 as it is a superset of GAF1
impl fmt::Display for Record {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> Result<(), fmt::Error> {
        let s = self
            .into_csv(&GafType::GAF2, ":".to_string(), "|".to_string())
            .join("\t");
        write!(f, "{}", s)
    }
}

impl Record {
    /// Create a new GAF record.
    pub fn new() -> Self {
        Record {
            db: String::new(),
            db_object_id: String::new(),
            db_object_symbol: String::new(),
            qualifier: Vec::<String>::new(),
            go_id: String::new(),
            db_ref: MultiMap::<String, String>::new(),
            evidence_code: String::new(),
            with_from: Vec::<String>::new(),
            aspect: String::new(),
            db_object_name: String::new(),
            db_object_synonym: Vec::<String>::new(),
            db_object_type: String::new(),
            taxon: Vec::<String>::new(),
            date: String::new(),
            assigned_by: String::new(),
            annotation_extension: None,
            gene_product_form_id: None,
        }
    }

    pub fn into_csv(
        &self,
        fileformat: &GafType,
        delimiter: String,
        terminator: String,
    ) -> Vec<String> {
        let db_ref = self
            .db_ref
            .iter()
            .map(|(a, b)| format!("{}{}{}", a, delimiter, b))
            .join(&terminator);

        let mut v = vec![
            self.db.clone(),
            self.db_object_id.clone(),
            self.db_object_symbol.clone(),
            self.qualifier.join(&terminator),
            self.go_id.clone(),
            db_ref,
            self.evidence_code.clone(),
            self.with_from.join(&terminator),
            self.aspect.clone(),
            self.db_object_name.clone(),
            self.db_object_synonym.join(&terminator),
            self.db_object_type.clone(),
            self.taxon.join(&terminator),
            self.date.clone(),
            self.assigned_by.clone(),
        ];

        if *fileformat == GafType::GAF2 {
            v.push(match &self.annotation_extension {
                Some(ae) => ae.join(&terminator),
                None => String::new(),
            });

            v.push(match &self.gene_product_form_id {
                Some(gene_product_form_id) => gene_product_form_id.clone(),
                None => String::new(),
            });
        }
        v
    }

    /// Get reference on db
    pub fn db(&self) -> &str {
        &self.db
    }

    /// Get reference on db_object_id
    pub fn db_object_id(&self) -> &str {
        &self.db_object_id
    }

    /// Get reference on db_object_symbol
    pub fn db_object_symbol(&self) -> &str {
        &self.db_object_symbol
    }

    /// Get reference on qualifier
    pub fn qualifier(&self) -> &Vec<String> {
        &self.qualifier
    }

    /// Get reference on go_id
    pub fn go_id(&self) -> &str {
        &self.go_id
    }

    /// Get reference on db_ref
    pub fn db_ref(&self) -> &MultiMap<String, String> {
        &self.db_ref
    }

    /// Get reference on evidence_code
    pub fn evidence_code(&self) -> &str {
        &self.evidence_code
    }

    /// Get reference on with_from
    pub fn with_from(&self) -> &Vec<String> {
        &self.with_from
    }

    /// Get reference on aspect
    pub fn aspect(&self) -> &str {
        &self.aspect
    }

    /// Get reference on db_object_name
    pub fn db_object_name(&self) -> &str {
        &self.db_object_name
    }

    /// Get reference on db_object_synonym
    pub fn db_object_synonym(&self) -> &Vec<String> {
        &self.db_object_synonym
    }

    /// Get reference on db_object_type
    pub fn db_object_type(&self) -> &str {
        &self.db_object_type
    }

    /// Get reference on taxon
    pub fn taxon(&self) -> &Vec<String> {
        &self.taxon
    }

    /// Get reference on date
    pub fn date(&self) -> &str {
        &self.date
    }

    /// Get reference on assigned_by
    pub fn assigned_by(&self) -> &str {
        &self.assigned_by
    }

    /// Get reference on annotation_extension
    pub fn annotation_extension(&self) -> &Option<Vec<String>> {
        &self.annotation_extension
    }

    /// Get reference on gene_product_form_id
    pub fn gene_product_form_id(&self) -> &Option<String> {
        &self.gene_product_form_id
    }

    /// Get mutable reference on db
    pub fn db_mut(&mut self) -> &mut String {
        &mut self.db
    }

    /// Get mutable reference on db_object_id
    pub fn db_object_id_mut(&mut self) -> &mut String {
        &mut self.db_object_id
    }

    /// Get mutable reference on db_object_symbol
    pub fn db_object_symbol_mut(&mut self) -> &mut String {
        &mut self.db_object_symbol
    }

    /// Get mutable reference on qualifier
    pub fn qualifier_mut(&mut self) -> &mut Vec<String> {
        &mut self.qualifier
    }

    /// Get mutable reference on go_id
    pub fn go_id_mut(&mut self) -> &mut String {
        &mut self.go_id
    }

    /// Get mutable reference on db_ref
    pub fn db_ref_mut(&mut self) -> &mut MultiMap<String, String> {
        &mut self.db_ref
    }

    /// Get mutable reference on evidence_code
    pub fn evidence_code_mut(&mut self) -> &mut String {
        &mut self.evidence_code
    }

    /// Get mutable reference on with_from
    pub fn with_from_mut(&mut self) -> &mut Vec<String> {
        &mut self.with_from
    }

    /// Get mutable reference on aspect
    pub fn aspect_mut(&mut self) -> &mut String {
        &mut self.aspect
    }

    /// Get mutable reference on db_object_name
    pub fn db_object_name_mut(&mut self) -> &mut String {
        &mut self.db_object_name
    }

    /// Get mutable reference on db_object_synonym
    pub fn db_object_synonym_mut(&mut self) -> &mut Vec<String> {
        &mut self.db_object_synonym
    }

    /// Get mutable reference on db_object_type
    pub fn db_object_type_mut(&mut self) -> &mut String {
        &mut self.db_object_type
    }

    /// Get mutable reference on taxon
    pub fn taxon_mut(&mut self) -> &mut Vec<String> {
        &mut self.taxon
    }

    /// Get mutable reference on date
    pub fn date_mut(&mut self) -> &mut String {
        &mut self.date
    }

    /// Get mutable reference on assigned_by
    pub fn assigned_by_mut(&mut self) -> &mut String {
        &mut self.assigned_by
    }

    /// Get mutable reference on annotation_extension
    pub fn annotation_extension_mut(&mut self) -> &mut Option<Vec<String>> {
        &mut self.annotation_extension
    }

    /// Get mutable reference on gene_product_form_id
    pub fn gene_product_form_id_mut(&mut self) -> &mut Option<String> {
        &mut self.gene_product_form_id
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const GAF_FILE1: &'static [u8] = b"!gaf-version: 1.0\n\
    GeneDB\tPF3D7_0100100.1\tVAR\t\tGO:0020002\tPMID:11420100\tTAS\t\tC\terythrocyte membrane protein 1, PfEMP1\tPFA0005w|VAR-UPSB1|MAL1P4.01\tmRNA\ttaxon:36329\t20190204\tGeneDB\n\
    GeneDB\tPF3D7_0100300.1\tVAR\t\tGO:0046789\tGO_REF:0000002\tIEA\tInterPro:IPR008602\tF\terythrocyte membrane protein 1, PfEMP1\tPFA0015c|VAR-UPSA3|MAL1P4.03|VAR3\tmRNA\ttaxon:36329\t20181008\tGeneDB\n";

    #[test]
    fn test_reader_gaf() {
        let data = vec![
            vec!["GeneDB", "GeneDB"],
            vec!["PF3D7_0100100.1", "PF3D7_0100300.1"],
            vec!["VAR", "VAR"],
            vec!["", ""],
            vec!["GO:0020002", "GO:0046789"],
            vec!["PMID:11420100", "GO_REF:0000002"],
            vec!["TAS", "IEA"],
            vec!["", "InterPro:IPR008602"],
            vec!["C", "F"],
            vec![
                "erythrocyte membrane protein 1, PfEMP1",
                "erythrocyte membrane protein 1, PfEMP1",
            ],
            vec![
                "PFA0005w|VAR-UPSB1|MAL1P4.01",
                "PFA0015c|VAR-UPSA3|MAL1P4.03|VAR3",
            ],
            vec!["mRNA", "mRNA"],
            vec!["taxon:36329", "taxon:36329"],
            vec!["20190204", "20181008"],
            vec!["GeneDB", "GeneDB"],
        ];

        let mut reader = Reader::new(GAF_FILE1, GafType::GAF1);
        for (i, r) in reader.records().enumerate() {
            let record = r.unwrap();

            let db_ref = record.db_ref();
            let mut db_ref_v: Vec<String> = Vec::<String>::new();
            for (key, values) in db_ref {
                for value in values {
                    let s = format!("{}:{}", key, value);
                    db_ref_v.push(s);
                }
            }

            assert_eq!(record.db(), data[0][i]);
            assert_eq!(record.db_object_id(), data[1][i]);
            assert_eq!(record.db_object_symbol(), data[2][i]);
            assert_eq!(record.qualifier().join("|"), data[3][i]);
            assert_eq!(record.go_id(), data[4][i]);
            assert_eq!(db_ref_v.join("|"), data[5][i]);
            assert_eq!(record.evidence_code(), data[6][i]);
            assert_eq!(record.with_from().join("|"), data[7][i]);
            assert_eq!(record.aspect(), data[8][i]);
            assert_eq!(record.db_object_name(), data[9][i]);
            assert_eq!(record.db_object_synonym().join("|"), data[10][i]);
            assert_eq!(record.db_object_type(), data[11][i]);
            assert_eq!(record.taxon().join("|"), data[12][i]);
            assert_eq!(record.date(), data[13][i]);
            assert_eq!(record.assigned_by(), data[14][i]);
            // No annotation_extension or gene_product_form_id, GFF1 doesn't have that!
        }
    }

    #[test]
    fn test_writer_gaf3() {
        let mut reader = Reader::new(GAF_FILE1, GafType::GAF1);
        let mut writer = Writer::new(vec![], GafType::GAF1);
        for r in reader.records() {
            writer
                .write(&r.ok().expect("Error reading record"))
                .ok()
                .expect("Error writing record");
        }

        let out = writer.inner.into_inner().unwrap();
        let s1 = String::from("!gaf-version: 1.0\n").to_owned()
            + String::from_utf8(out).unwrap().as_str();
        let s2 = String::from_utf8(GAF_FILE1.to_owned()).unwrap();
        assert_eq!(s1, s2)
    }

    #[test]
    fn test_gaftype_separator() {
        assert_eq!(GafType::GAF1.separator(), (b':', b'|', 0u8));
        assert_eq!(GafType::GAF2.separator(), (b':', b'|', 0u8));
    }

    #[test]
    fn test_gaftype_name() {
        assert_eq!(GafType::GAF1.name(), "GAF1");
        assert_eq!(GafType::GAF2.name(), "GAF2");
    }

    #[test]
    fn test_gaf_record_db() {
        let mut r = Record::new();
        assert_eq!(r.db(), "");
        *r.db_mut() = "TEST".to_string();
        assert_eq!(r.db(), "TEST");
    }

    #[test]
    fn test_gaf_record_db_object_id() {
        let mut r = Record::new();
        assert_eq!(r.db_object_id(), "");
        *r.db_object_id_mut() = "TEST".to_string();
        assert_eq!(r.db_object_id(), "TEST");
    }

    #[test]
    fn test_gaf_record_db_object_symbol() {
        let mut r = Record::new();
        assert_eq!(r.db_object_symbol(), "");
        *r.db_object_symbol_mut() = "TEST".to_string();
        assert_eq!(r.db_object_symbol(), "TEST");
    }

    #[test]
    fn test_gaf_record_go_id() {
        let mut r = Record::new();
        assert_eq!(r.go_id(), "");
        *r.go_id_mut() = "TEST".to_string();
        assert_eq!(r.go_id(), "TEST");
    }

    #[test]
    fn test_gaf_record_evidence_code() {
        let mut r = Record::new();
        assert_eq!(r.evidence_code(), "");
        *r.evidence_code_mut() = "TEST".to_string();
        assert_eq!(r.evidence_code(), "TEST");
    }

    #[test]
    fn test_gaf_record_aspect() {
        let mut r = Record::new();
        assert_eq!(r.aspect(), "");
        *r.aspect_mut() = "TEST".to_string();
        assert_eq!(r.aspect(), "TEST");
    }

    #[test]
    fn test_gaf_record_date() {
        let mut r = Record::new();
        assert_eq!(r.date(), "");
        *r.date_mut() = "TEST".to_string();
        assert_eq!(r.date(), "TEST");
    }

    #[test]
    fn test_gaf_record_assigned_by() {
        let mut r = Record::new();
        assert_eq!(r.assigned_by(), "");
        *r.assigned_by_mut() = "TEST".to_string();
        assert_eq!(r.assigned_by(), "TEST");
    }

    #[test]
    fn test_gaf_record_db_object_name() {
        let mut r = Record::new();
        assert_eq!(r.db_object_name(), "");
        *r.db_object_name_mut() = "TEST".to_string();
        assert_eq!(r.db_object_name(), "TEST");
    }

    #[test]
    fn test_gaf_record_db_object_type() {
        let mut r = Record::new();
        assert_eq!(r.db_object_type(), "");
        *r.db_object_type_mut() = "TEST".to_string();
        assert_eq!(r.db_object_type(), "TEST");
    }

    #[test]
    fn test_gaf_record_qualifier() {
        let mut r = Record::new();
        let test = vec!["TEST".to_string()];
        assert!(r.qualifier().is_empty());
        *r.qualifier_mut() = test.clone();
        assert_eq!(*r.qualifier(), test);
    }

    #[test]
    fn test_gaf_record_with_from() {
        let mut r = Record::new();
        let test = vec!["TEST".to_string()];
        assert!(r.with_from().is_empty());
        *r.with_from_mut() = test.clone();
        assert_eq!(*r.with_from(), test);
    }

    #[test]
    fn test_gaf_record_db_object_synonym() {
        let mut r = Record::new();
        let test = vec!["TEST".to_string()];
        assert!(r.db_object_synonym().is_empty());
        *r.db_object_synonym_mut() = test.clone();
        assert_eq!(*r.db_object_synonym(), test);
    }

    #[test]
    fn test_gaf_record_taxon() {
        let mut r = Record::new();
        let test = vec!["TEST".to_string()];
        assert!(r.taxon().is_empty());
        *r.taxon_mut() = test.clone();
        assert_eq!(*r.taxon(), test);
    }

    #[test]
    fn test_gaf_record_annotation_extension() {
        let mut r = Record::new();
        let test = Some(vec!["TEST".to_string()]);
        assert!(r.annotation_extension().is_none());
        *r.annotation_extension_mut() = test.clone();
        assert_eq!(*r.annotation_extension(), test);
    }

    #[test]
    fn test_gaf_record_gene_product_form_id() {
        let mut r = Record::new();
        let test = Some("TEST".to_string());
        assert!(r.gene_product_form_id().is_none());
        *r.gene_product_form_id_mut() = test.clone();
        assert_eq!(*r.gene_product_form_id(), test);
    }

    #[test]
    fn test_gaf_record_db_ref() {
        let mut r = Record::new();
        let mut test = MultiMap::new();
        test.insert("k1".to_string(), "v1".to_string());
        test.insert("k1".to_string(), "v2".to_string());
        test.insert("k3".to_string(), "v3".to_string());
        assert!(r.db_ref().is_empty());
        *r.db_ref_mut() = test.clone();
        assert_eq!(*r.db_ref(), test);
    }

    #[test]
    fn test_gaf_record_fmt() {
        let line = "GeneDB\tPF3D7_0100100.1\tVAR\t\tGO:0020002\tPMID:11420100\tTAS\t\tC\terythrocyte membrane protein 1, PfEMP1\tPFA0005w|VAR-UPSB1|MAL1P4.01\tmRNA\ttaxon:36329\t20190204\tGeneDB\t\t";
        let mut reader = Reader::new(GAF_FILE1, GafType::GAF1);
        for (_i, r) in reader.records().enumerate() {
            let record = r.unwrap();
            assert_eq!(format!("{}", record), line);
            break;
        }
    }
}

#![allow(dead_code)]

use anyhow::Error;
// use left_align_indel::left_align_indels;
use minimap2::{Aligner, Built};
use rust_htslib::bam::{Header, HeaderView, Record};

use crate::alignment::AlignmentRef;

pub type Remapper = Aligner<Built>;

const DUMMY_REFERENCE: &[u8; 3] = b"ACT";

pub struct Realigner {
    aligner: Remapper,
    refname: Vec<u8>,
    headerview: Option<HeaderView>,
}

/// remove supplementary maps, get only one alignment per read.
/// this is meant to be called on a [Vec<Mapping>] from a single read
pub fn filter_maps(maps: &mut Vec<Record>) -> Record {
    let mut aln: Vec<Record>;

    aln = maps
        .drain(..)
        .filter_map(|m| {
            if !m.is_supplementary() && !m.is_secondary() {
                Some(m)
            } else {
                None
            }
        })
        .collect();

    assert_eq!(aln.len(), 1);
    aln.remove(0)
}

pub enum AlignerReference<'a> {
    Sequence(&'a [u8]),
    Fasta(&'a str),
}

impl Realigner {
    pub fn init_to_ref(
        &mut self,
        reference: AlignerReference,
        refname: Option<&str>,
    ) -> Result<(), Error> {
        let aligner: Aligner<Built>;
        let mut aligner_build = Aligner::builder()
            .sr()
            .with_cigar()
            .with_sam_out()
            .with_index_threads(num_cpus::get());

        aligner_build.mapopt.best_n = 1;
        aligner_build.mapopt.b = 6;
        aligner_build.mapopt.q = 12;
        aligner_build.mapopt.end_bonus = 100;
        aligner_build.idxopt.k = 5;
        aligner_build.idxopt.w = 5;

        match reference {
            AlignerReference::Fasta(file) => {
                aligner = aligner_build.with_index(file, None).map_err(Error::msg)?
            }
            AlignerReference::Sequence(bytes) => {
                aligner = aligner_build.with_seq(bytes).map_err(Error::msg)?
            }
        };

        let refname = refname.unwrap_or("REF").as_bytes().to_vec();

        let mut header = Header::new();
        aligner.populate_header(&mut header);
        let headerview = HeaderView::from_header(&header);

        self.aligner = aligner;
        self.refname = refname;
        self.headerview = Some(headerview);

        Ok(())
    }

    pub fn build_empty() -> Result<Self, Error> {
        let aligner = Aligner::builder().with_seq(b"ACT").map_err(Error::msg)?;

        Ok(Self {
            aligner,
            refname: "NONE".into(),
            headerview: None,
        })
    }

    pub fn realign(
        &self,
        rec: &Record,
        header: &HeaderView,
        outvec: &mut Vec<Record>,
    ) -> Result<(), Error> {
        *outvec = self
            .aligner
            .map_to_sam(
                rec.seq().as_bytes().as_slice(),
                None,
                Some(rec.qname()),
                header,
                None,
                None,
            )
            .expect("failed to realign!");

        Ok(())
    }

    pub fn realign_region_record(&mut self, records: &mut [Record]) -> Result<(), Error> {
        let mut aln: Record;
        let mut maps: Vec<Record> = vec![];
        let header = self.headerview.as_ref().unwrap();

        for rec in records.iter_mut() {
            self.realign(rec, header, &mut maps)?;
            if maps.is_empty() {
                continue;
            }

            aln = filter_maps(&mut maps);
            *rec = aln;
        }

        Ok(())
    }

    pub fn realign_region_plp(&mut self, pileups: &mut [AlignmentRef]) -> Result<(), Error> {
        let mut aln: Record;
        let mut maps: Vec<Record> = vec![];
        let header = self.headerview.as_ref().unwrap();

        for plp in pileups.iter_mut() {
            let plp = &mut plp.borrow_mut();
            self.realign(&plp.rec, header, &mut maps)?;
            if maps.is_empty() {
                continue;
            }

            aln = filter_maps(&mut maps);
            // record = aln;
            //
            plp.cstate.cig = aln.cigar().clone();
            plp.cstate.bam_pos = aln.pos() as u32;
            plp.rec = aln;
            plp.cstate.icig = 0;
            plp.cstate.iseq = 0;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;

    use std::fs::File;
    use std::io::{BufRead, BufReader};

    use rust_htslib::bam::Record;

    fn get_test_dir() -> PathBuf {
        let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
        path.push("test_data");
        path
    }

    fn get_test_file(test_file: &str) -> PathBuf {
        let mut path = get_test_dir();
        path.push(test_file);
        path
    }

    fn get_fastq_sequences(test_file: &str) -> Vec<Vec<u8>> {
        let mut out = vec![];
        let f = get_test_file(test_file);
        let mut lines = BufReader::new(File::open(f).unwrap()).lines();

        while let (Some(_header), Some(seq), Some(_plus), Some(_qual)) =
            (lines.next(), lines.next(), lines.next(), lines.next())
        {
            out.push(seq.unwrap().as_bytes().into());
        }

        out
    }

    fn bam_from_fastq(test_file: &str) -> Vec<Record> {
        let reads = get_fastq_sequences(test_file);
        let records = reads
            .iter()
            .map(|r| {
                let mut record = Record::new();
                record.set(b"4", None, r, vec![255_u8; r.len()].as_slice());
                record
            })
            .collect();

        records
    }

    #[test]
    fn test1() {
        let ref_file = get_test_file("cDNA.fasta").to_str().unwrap().to_string();
        let mut ref_records = bam_from_fastq("cDNA_reads.fq");

        let mut realigner = Realigner::build_empty().unwrap();

        realigner
            .init_to_ref(AlignerReference::Fasta(&ref_file), None)
            .unwrap();

        realigner.realign_region_record(&mut ref_records).unwrap();

        println! {"{:?}", ref_records}
        for r in &ref_records {
            assert!(!r.is_unmapped());
            println! {"{:?}, {:?}", r, r.cigar()}
        }
    }

    #[test]
    fn test2() {
        let ref_file = get_test_file("hiv.fasta").to_str().unwrap().to_string();
        let mut ref_records = bam_from_fastq("hiv_reads.fq");

        let mut realigner = Realigner::build_empty().unwrap();

        realigner
            .init_to_ref(AlignerReference::Fasta(&ref_file), None)
            .unwrap();

        realigner.realign_region_record(&mut ref_records).unwrap();

        println! {"{:?}", ref_records}
        for r in &ref_records {
            println! {"{:?}, {:?}", r, r.cigar()}
            assert!(!r.is_unmapped());
        }
    }
}

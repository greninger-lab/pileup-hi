#![allow(dead_code)]

use bio::alignment::{pairwise::Aligner, Alignment};
use bio::scores::blosum62;
use rust_htslib::bam::{
    record::{Cigar, CigarString},
    Record,
};

use crate::pileup_iterator::UNINIT_POS;

pub type Remapper = Aligner<fn(u8, u8) -> i32>;

pub struct Realigner {
    prev_w_start: i64,
    prev_w_end: i64,
    aligner: Remapper,
}

fn parse_cigar_string(cigar_str: &str) -> CigarString {
    let mut ops = Vec::new();
    let mut current_num = 0;
    for c in cigar_str.chars() {
        match c {
            '0'..='9' => current_num = current_num * 10 + (c as u32 - '0' as u32),
            _ => {
                let op = match c {
                    'M' => Cigar::Match(current_num),
                    'I' => Cigar::Ins(current_num),
                    'D' => Cigar::Del(current_num),
                    'N' => Cigar::RefSkip(current_num),
                    'S' => Cigar::SoftClip(current_num),
                    'H' => Cigar::HardClip(current_num),
                    'P' => Cigar::Pad(current_num),
                    'X' => Cigar::Diff(current_num),
                    '=' => Cigar::Equal(current_num),
                    _ => panic!("Invalid CIGAR operation: {}", c),
                };
                ops.push(op);
                current_num = 0;
            }
        }
    }

    CigarString(ops)
}

pub struct RefWindow<'a> {
    seq: &'a [u8],
    start_offset: i64,
}

impl Realigner {
    pub fn new() -> Self {
        let aligner: Remapper = Aligner::new(-5, -1, blosum62);

        Self {
            prev_w_start: 0,
            prev_w_end: UNINIT_POS,
            aligner,
        }
    }

    /// update the min and max window to avoid repeating work.
    pub fn realign_region(&mut self, ref_window: RefWindow, mut records: Vec<&mut Record>) {
        let mut aln: Alignment = Alignment::default();

        records.iter_mut().for_each(|rec| {
            aln = self
                .aligner
                .semiglobal(ref_window.seq, &rec.seq().as_bytes());
            rec.set_cigar(Some(&parse_cigar_string(&aln.cigar(false))));
            rec.set_pos(aln.xstart as i64 + ref_window.start_offset);
            println!{"{}\n YSTART: {}\n XSTART: {}", aln.pretty(ref_window.seq, &rec.seq().as_bytes(), 50), aln.ystart, aln.xstart}
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::Record;

    #[test]
    fn test1() {
        let mut record = Record::new();
        record.set(
            b"read1", None, b"ATCATT", b"######",
            //012345
        );

        record.set_pos(10);

        let seq = b"GATCATTATGATAT";
        let rw = RefWindow {
            seq: &seq[2..],
            start_offset: 2,
        };

        let mut realigner = Realigner::new();
        realigner.realign_region(rw, vec![&mut record]);
        assert_eq!(record.pos(), 2);
    }
}

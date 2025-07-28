use rust_htslib::bam::{
    record::{Cigar, CigarStringView},
    Record,
};

pub struct CigarState {
    pub cig: CigarStringView,
    pub icig: usize,  // position in cigar string
    pub iseq: u32,    // position in read sequence that corresponds to cigar pos
    pub bam_pos: u32, // ref coord of first base
    pub qpos: usize,
    pub del: bool,
}

pub struct PileUp {
    pub rec: Record,
    pub cstate: CigarState,
}

pub fn cigar2rlen(r: &Record) -> i64 {
    let mut len = 0;
    for op in &r.cigar() {
        match op {
            Cigar::Match(l)
            | Cigar::Del(l)
            | Cigar::RefSkip(l)
            | Cigar::Equal(l)
            | Cigar::Diff(l) => len += l,
            _ => (),
        }
    }

    len as i64
}

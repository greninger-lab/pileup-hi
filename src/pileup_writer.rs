use crate::pileup_iterator::CigarState;
use crate::refseq::RefSeq;

use anyhow::Error;
use rust_htslib::bam::{ext::BamRecordExtensions, record::Cigar, Record};
use std::io::Write;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

pub fn get_base(mut cur_base: u8, is_reverse: bool) -> u8 {
    match is_reverse {
        false => cur_base.make_ascii_uppercase(),
        true => cur_base.make_ascii_lowercase(),
    }

    cur_base
}

// cap qualitites at max of 126; this also helps avoid non-ascii output
pub fn get_qual(qual: u8) -> u8 {
    match qual.cmp(&92).is_gt() {
        true => 126,
        false => qual + 33,
    }
}

pub fn get_base_to_ref(
    cur_base: u8,
    plp_ref_pos: i64,
    refseq: &Option<RefSeq>,
    is_reverse: bool,
) -> Result<u8, Error> {
    if let Some(refseq) = refseq {
        if refseq.get_base(plp_ref_pos as u64)? == cur_base {
            if is_reverse {
                return Ok(R_MATCH);
            } else {
                return Ok(F_MATCH);
            }
        }
    }
    return Ok(get_base(cur_base, is_reverse));
}

pub struct PileupWriter {
    qual_buf: Vec<u8>,
    seq_buf: Vec<u8>,
}

impl PileupWriter {
    pub fn write_match(
        &mut self,
        rec: &Record,
        plp_read_idx: u32,
        plp_ref_pos: i64,
        refseq: &Option<RefSeq>,
    ) -> Result<(), Error> {
        if plp_ref_pos == rec.pos() {
            self.seq_buf.push(FIRST_POS);
            self.seq_buf.push(get_qual(rec.mapq()));
        }

        let base = get_base_to_ref(
            rec.seq()[plp_read_idx as usize],
            plp_ref_pos,
            refseq,
            rec.is_reverse(),
        )?;
        let qual = get_qual(rec.qual()[plp_read_idx as usize]);

        self.seq_buf.push(base);
        self.qual_buf.push(qual);

        if plp_ref_pos == rec.reference_end() - 1 {
            self.seq_buf.push(LAST_POS);
        }

        Ok(())
    }

    pub fn write_deletion(&mut self, qual: u8) {
        self.seq_buf.push(b'*');
        self.qual_buf.push(get_qual(qual));
    }

    pub fn write_insertion(
        &mut self,
        cs: &CigarState,
        r: &Record,
        plp_read_pos: u32,
    ) -> Result<(), Error> {
        let mut k = cs.icig + 1; // move into insertion
        let ncig = cs.cig.len();
        while k < ncig {
            match cs.cig[k] {
                Cigar::Pad(l) => {
                    self.seq_buf.extend(std::iter::repeat_n(b'*', l as usize));
                }

                Cigar::Ins(l) => {
                    write!(self.seq_buf, "+{}", l)?;

                    // starting at plp_read_pos +1 so we begin at the first base of the insertion
                    let (s, e) = ((plp_read_pos + 1) as usize, (plp_read_pos + 1 + l) as usize);
                    for i in s..e {
                        let base = get_base(r.seq()[i], r.is_reverse());
                        self.seq_buf.push(base);
                    }
                }

                _ => break,
            }

            k += 1;
        }

        Ok(())
    }

    pub fn write_deletion_start(
        &mut self,
        rec: &Record,
        plp_ref_pos: i64,
        refseq: &Option<RefSeq>,
        del_len: i64,
    ) -> Result<(), Error> {
        write!(self.seq_buf, "-{}", del_len)?;

        for p in plp_ref_pos..plp_ref_pos + del_len {
            let b = match refseq {
                Some(refseq) => get_base(refseq.get_base(p as u64)?, rec.is_reverse()),
                None => b'N',
            };

            self.seq_buf.push(get_base(b, rec.is_reverse()));
        }

        Ok(())
    }

    pub fn write_pileup_str(
        &mut self,
        ref_base: u8,
        plp_ref_pos: i64,
        nbases: usize,
        nins: usize,
        ndel: usize,
        tidname: &str,
    ) -> Result<(), Error> {
        print! {"{}\t{}\t{}\t{}\t", tidname, plp_ref_pos + 1, char::from(ref_base), nbases + nins + ndel }

        if self.seq_buf.is_empty() {
            print! {"*\t"}
        } else {
            print! {"{}\t", std::str::from_utf8(&self.seq_buf)?}
            self.seq_buf.clear();
        }

        if self.qual_buf.is_empty() {
            print! {"*"}
        } else {
            print! {"{}", std::str::from_utf8(&self.qual_buf)?}
            self.qual_buf.clear();
        }

        print! {"\n"}

        Ok(())
    }

    pub fn new() -> Self {
        let (seq_buf, qual_buf) = (Vec::with_capacity(500), Vec::with_capacity(500));

        Self { seq_buf, qual_buf }
    }
}

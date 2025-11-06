use crate::{
    alignment::CigarState,
    output::OrderedPileupOutput,
    pileup_iterator::{PileupBaseCall, PileupPayload},
};
use crossbeam::channel::Sender;

use anyhow::Error;
use rust_htslib::bam::{ext::BamRecordExtensions, record::Cigar, Record};
use std::io::Write;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

const F_REFSKIP: u8 = b'>';
const R_REFSKIP: u8 = b'<';

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

pub fn get_base_to_ref(cur_base: u8, ref_base: u8, is_reverse: bool) -> Result<u8, Error> {
    if ref_base == cur_base {
        if is_reverse {
            return Ok(R_MATCH);
        } else {
            return Ok(F_MATCH);
        }
    }
    Ok(get_base(cur_base, is_reverse))
}

pub type PileupStringOutput<T: OrderedPileupOutput> = Sender<T>;

pub struct PileupWriter {
    inner: PileupStringType,
}

impl PileupWriter {
    pub fn new_inplace() -> Self {
        Self {
            inner: PileupStringType::InPlace(PileupStringInPlace {
                plp_string: PileupString::new(),
            }),
        }
    }

    pub fn new_multi(s: Sender<PileupString>) -> Self {
        Self {
            inner: PileupStringType::MultiThreaded(PileupStringMultiThread {
                plp_string: PileupString::new(),
                out: s,
            }),
        }
    }

    pub fn intake(&mut self, p: PileupPayload) -> Result<(), Error> {
        match &mut self.inner {
            PileupStringType::InPlace(ref mut s) => s.plp_string.intake(p),
            PileupStringType::MultiThreaded(ref mut s) => s.plp_string.intake(p),
        }
    }

    pub fn write(&mut self) -> Result<(), Error> {
        match &mut self.inner {
            PileupStringType::InPlace(ref mut s) => s.plp_string.write_pileup_str(),
            PileupStringType::MultiThreaded(ref mut s) => {
                s.out.send(s.plp_string.clone()).map_err(Error::msg)?;
                s.plp_string.qual_buf.shrink_to(0);
                s.plp_string.seq_buf.shrink_to(0);
                Ok(())
            }
        }
    }

    pub fn update(&mut self, tid: i32, ref_pos: i64, ref_base: u8, ref_name: String, depth: u32) {
        match &mut self.inner {
            PileupStringType::InPlace(ref mut s) => {
                s.plp_string.update(tid, ref_pos, ref_base, ref_name, depth)
            }
            PileupStringType::MultiThreaded(ref mut s) => {
                s.plp_string.update(tid, ref_pos, ref_base, ref_name, depth)
            }
        }
    }
}

/// Class with methods to write pileup information output to stdout/file that is meant to be
/// compliant with the default format of samtools mpileup.
#[derive(Clone)]
pub struct PileupString {
    tid: i32,
    ref_pos: i64,
    ref_base: u8,
    ref_name: String,
    depth: u32,
    seq_buf: Vec<u8>,
    qual_buf: Vec<u8>,
}

pub enum PileupStringType {
    InPlace(PileupStringInPlace),
    MultiThreaded(PileupStringMultiThread),
}

/// A pileup string that is not copied, and simply printed to stdout once it is considered full.
/// Should not be used when multiple strings are being developed concurrently.
pub struct PileupStringInPlace {
    pub plp_string: PileupString,
}

pub struct PileupStringMultiThread {
    pub plp_string: PileupString,
    out: Sender<PileupString>,
}

impl PileupString {
    pub fn intake(&mut self, p: PileupPayload) -> Result<(), Error> {
        match p.call {
            PileupBaseCall::Match => {
                self.write_match(p.rec, p.plp_read_pos as u32, p.plp_ref_pos, p.ref_base)?;
            }

            PileupBaseCall::DeletionStart => {
                self.write_match(p.rec, p.plp_read_pos as u32, p.plp_ref_pos, p.ref_base)?;
                self.write_deletion_start(p.rec, p.aux, p.aux.len() as i64)?;
            }

            PileupBaseCall::RefSkip => {
                self.write_refskip(p.rec, p.plp_read_pos as u32, p.plp_ref_pos, p.ref_base)?;
            }

            PileupBaseCall::Insertion => {
                self.write_match(p.rec, p.plp_read_pos as u32, p.plp_ref_pos, p.ref_base)?;
                self.write_insertion(p.cstate, p.rec, p.plp_read_pos as u32)?;
            }

            PileupBaseCall::Gap => {
                self.write_deletion(*p.rec.qual().get(p.plp_read_pos).unwrap_or(&0))
            }

            PileupBaseCall::NA => (),
        }

        Ok(())
    }

    pub fn write_refskip(
        &mut self,
        rec: &Record,
        plp_read_idx: u32,
        _plp_ref_pos: i64,
        _ref_base: u8,
    ) -> Result<(), Error> {
        let qual = get_qual(rec.qual()[plp_read_idx as usize]);

        let sym = if rec.is_reverse() {
            R_REFSKIP
        } else {
            F_REFSKIP
        };

        // self.seq_buf.push(base);
        self.seq_buf.push(sym);
        self.qual_buf.push(qual);

        Ok(())
    }

    pub fn write_match(
        &mut self,
        rec: &Record,
        plp_read_idx: u32,
        plp_ref_pos: i64,
        ref_base: u8,
    ) -> Result<(), Error> {
        if plp_ref_pos == rec.pos() {
            self.seq_buf.push(FIRST_POS);
            self.seq_buf.push(get_qual(rec.mapq()));
        }

        let base = get_base_to_ref(rec.seq()[plp_read_idx as usize], ref_base, rec.is_reverse())?;
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
        del_seq: &[u8],
        del_len: i64,
    ) -> Result<(), Error> {
        write!(self.seq_buf, "-{}", del_len)?;
        del_seq
            .iter()
            .for_each(|b| self.seq_buf.push(get_base(*b, rec.is_reverse())));
        Ok(())
    }

    pub fn write_pileup_str(&mut self) -> Result<(), Error> {
        print! {"{}\t{}\t{}\t{}\t", self.ref_name, self.ref_pos + 1, char::from(self.ref_base), self.depth }

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

        println! {};

        Ok(())
    }

    pub fn new() -> Self {
        let (seq_buf, qual_buf) = (Vec::with_capacity(500), Vec::with_capacity(500));

        Self {
            tid: 0,
            ref_pos: 0,
            ref_base: 0,
            ref_name: "".to_string(),
            depth: 0,
            seq_buf,
            qual_buf,
        }
    }

    pub fn update(&mut self, tid: i32, ref_pos: i64, ref_base: u8, ref_name: String, depth: u32) {
        self.tid = tid;
        self.ref_pos = ref_pos;
        self.ref_base = ref_base;
        self.ref_name = ref_name;
        self.depth = depth;
    }
}

impl OrderedPileupOutput for PileupString {
    fn tid(&self) -> i32 {
        self.tid
    }

    fn pos(&self) -> i64 {
        self.ref_pos
    }

    fn write(&mut self) -> Result<(), Error> {
        self.write_pileup_str()
    }
}

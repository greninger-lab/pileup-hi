use crate::alignment::PileupAlignment;
use crate::errors::{Error, ErrorKind};
use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::{cell::RefCell, rc::Rc};

extern "C" {
    fn tweak_overlap_quality(a: *mut rust_htslib::htslib::bam1_t, b: *mut rust_htslib::htslib::bam1_t) -> i32;
}

pub fn tweak_overlap_qual(a: &mut Record, b: &mut Record) -> Result<(), Error> {
    unsafe {
        let ret = tweak_overlap_quality(a.inner_mut() as *mut _, b.inner_mut() as *mut _);
        if ret < 0 {
            let qname = std::str::from_utf8(a.qname())?;
            return Err(Error::from(ErrorKind::MateOverlapFailed(qname.to_string())));
        }
    }
    Ok(())
}

pub type OverlapMap = HashMap<u64, Rc<RefCell<PileupAlignment>>>;

pub trait MapOverlaps {
    fn push(&mut self, r: Rc<RefCell<PileupAlignment>>);
    fn delete_hash(&mut self, r: u64);
    fn delete_read(&mut self, r: &Record);
}

pub fn hash_qname(r: &Record) -> u64 {
    let mut hasher = DefaultHasher::new();
    r.qname().hash(&mut hasher);
    hasher.finish()
}

impl MapOverlaps for OverlapMap {
    fn push(&mut self, plp: Rc<RefCell<PileupAlignment>>) {
        let mut _r = plp.borrow_mut();
        let len = _r.cstate.read_len_from_cigar;
        let r = &mut _r.rec;

        if r.is_mate_unmapped() || !r.is_proper_pair() {
            return;
        }

        if (r.mtid() >= 0 && (r.mtid() != r.tid()))
            || r.insert_size().abs() >= 2 * (r.seq_len() as i64) && r.mpos() >= r.pos() + len
        {
            return;
        }

        let h = hash_qname(r);

        if let Some(mate) = self.get_mut(&h) {
            tweak_overlap_qual(&mut mate.borrow_mut().rec, r).unwrap();
            self.delete_hash(h);
        } else if r.mpos() >= r.pos() || (r.is_paired() && r.mpos() == -1) {
            // criteria passed, insert
            self.insert(h, Rc::clone(&plp));
        }
    }

    fn delete_hash(&mut self, r: u64) {
        self.remove(&r);
    }

    fn delete_read(&mut self, r: &Record) {
        let mut h = DefaultHasher::new();
        r.qname().hash(&mut h);
        self.remove(&h.finish());
    }
}

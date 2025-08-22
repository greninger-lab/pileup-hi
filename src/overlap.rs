#![allow(dead_code)]
// use crate::alignment::Pileup;
use crate::alignment::Alignment;
use crate::read_walker::WalkMatches;
use anyhow::Error;
use rust_htslib::bam::Record;
use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::slice;
use std::{cell::RefCell, cmp::Ordering, rc::Rc};

pub type OverlapMap = HashMap<u64, Rc<RefCell<Alignment>>>;

pub trait MapOverlaps {
    fn push(&mut self, r: Rc<RefCell<Alignment>>);
    fn delete_hash(&mut self, r: u64);
    fn delete_read(&mut self, r: &Record);
}

pub fn hash_qname(r: &Record) -> u64 {
    let mut hasher = DefaultHasher::new();
    r.qname().hash(&mut hasher);
    hasher.finish()
}

impl MapOverlaps for OverlapMap {
    fn push(&mut self, plp: Rc<RefCell<Alignment>>) {
        let r = &mut plp.borrow_mut().rec;

        if r.is_mate_unmapped() || !r.is_proper_pair() {
            return;
        }

        if r.mtid() >= 0 && (r.mtid() != r.tid()) {
            return;
        }

        let h = hash_qname(r);

        if let Some(mate) = self.get_mut(&h) {
            tweak_overlap_qual(&mut mate.borrow_mut().rec, r).unwrap();
            self.delete_hash(h);
            return;
        }

        if r.mpos() >= r.pos() || (r.is_paired() && r.mpos() == -1) {
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

/// this is just a modified version of [Record::qual] that returns a mutable slice to the qual
/// array.
///
/// check to make sure div_ceil is ok
pub fn qual_mut(r: &mut Record) -> &mut [u8] {
    let i: usize = r.inner.core.l_qname as usize + r.cigar_len() * 4 + r.seq_len().div_ceil(2);
    let dat_mut = unsafe { slice::from_raw_parts_mut(r.inner().data, r.inner().l_data as usize) };
    &mut dat_mut[i..][..r.seq_len()]
}

/// modify the phred score at a given position with a given value
pub fn set_qual(r: &mut Record, idx: usize, qual: u8) -> Result<(), Error> {
    if idx >= r.seq_len() {
        anyhow::bail!(
            "qual index of {} is out of bounds for record of len {}",
            idx,
            r.seq_len()
        )
    } else {
        qual_mut(r)[idx] = qual;
        Ok(())
    }
}

// designed to mimic arbitrary pseduorandom read selection from htslib's Wang hashmap
// we just combine the two rounds of hashing and return true/false based on first bit
// of hashed value.
pub fn decide_which_read(chars: &[u8]) -> bool {
    let mut h = chars[0] as u32;

    if h > 0 {
        for c in &chars[1..] {
            h = (h << 5) - (h) + *c as u32;
        }
    }

    h += !(h << 15);
    h ^= h >> 10;
    h += h << 3;
    h ^= h >> 6;
    h += !(h << 11);
    h ^= h >> 16;
    h & 1 != 0
}

/// at a reference position covered by both reads, null one's phred scores and keep the other
pub fn null_ref_bases(
    a: &mut Record,
    aread: usize,
    b: &mut Record,
    bread: usize,
    amul: bool,
    base_a: &mut u8,
    base_b: &mut u8,
    new_qual: &mut u8,
) -> Result<(), Error> {
    (*base_a, *base_b) = (a.seq()[aread], b.seq()[bread]);

    // both bases mismatch, so pick based on quality (or random if tie)
    if base_a != base_b {
        match a.qual()[aread].cmp(&b.qual()[bread]) {
            Ordering::Less => {
                *new_qual = (b.qual()[bread] as f64 * 0.8) as u8;
                set_qual(a, aread, 0)?;
                set_qual(b, bread, *new_qual)?;
            }

            Ordering::Greater => {
                *new_qual = (a.qual()[aread] as f64 * 0.8) as u8;
                set_qual(a, aread, *new_qual)?;
                set_qual(b, bread, 0)?;
            }

            Ordering::Equal => {
                if amul {
                    *new_qual = (a.qual()[aread] as f64 * 0.8) as u8;
                    set_qual(a, aread, *new_qual)?;
                    set_qual(b, bread, 0)?;
                } else {
                    *new_qual = (b.qual()[bread] as f64 * 0.8) as u8;
                    set_qual(a, aread, 0)?;
                    set_qual(b, bread, *new_qual)?;
                }
            }
        }

        // println! {"Adjusting mismatch: {aread} {bread} {} {} {}", a.qual()[aread], b.qual()[bread], std::str::from_utf8(a.qname())?}
    } else {
        // both bases match; Bump the quality up for one read and null the
        // other's
        *new_qual = (a.qual()[aread].wrapping_add(b.qual()[bread])).min(200);

        // set quals accordingly
        if amul {
            set_qual(a, aread, *new_qual)?;
            set_qual(b, bread, 0)?;
        } else {
            set_qual(a, aread, 0)?;
            set_qual(b, bread, *new_qual)?;
        }

        // println!("Adjusting quality to be ultra-confident A POS: {aread} | B POS: {bread} | A QUAL: {} | B QUAL: {} {}", a.qual()[aread], b.qual()[bread], std::str::from_utf8(a.qname())?)
    }

    Ok(())
}

pub fn tweak_overlap_qual(a: &mut Record, b: &mut Record) -> Result<(), Error> {
    let mut new_qual: u8 = 0;
    let mut base_a @ mut base_b: u8 = b'N';
    let amul @ bmul: bool;

    if a.pos() > b.pos() {
        std::mem::swap(a, b);
    }

    // we assume that we encounter reads in order (e.g coord-sorted).
    assert!(a.pos() <= b.pos());

    let mut ap = a.walk_matches();
    let mut bp = b.walk_matches();

    // catch up read A to read B's start pos
    // if no overlap, return early.

    // for now: using htslib's hashing heuristic to decide which read to modify, bound to change
    // (maybe)
    match decide_which_read(a.qname()) {
        true => (amul, bmul) = (true, false),
        false => (amul, bmul) = (false, true),
    }

    let mut iref = b.pos();
    let (mut apos, mut a_iref, mut bpos, mut b_iref) =
        (ap.read_pos, ap.genome_pos, bp.read_pos, bp.genome_pos);

    loop {
        while ap.genome_pos < iref {
            match ap.next() {
                Some((ap, ai)) => (apos, a_iref) = (ap, ai),
                None => return Ok(()),
            }
        }

        while bp.genome_pos < iref {
            match bp.next() {
                Some((bp, bi)) => (bpos, b_iref) = (bp, bi),
                None => return Ok(()),
            }
        }

        // print!("{iref}");

        iref = iref.max(ap.genome_pos);
        iref = iref.max(bp.genome_pos);
        iref += 1; // prepare for next position

        // check for deletion in read A
        // println! {"APOS: {apos} BPOS: {bpos} {} {} {} {}", ap.after_del(), bp.after_del(), ap.genome_pos, bp.genome_pos} // if a_iref > b_iref && ap.passed_deletion() {
        if a_iref > b_iref && ap.after_del() {
            while b_iref < a_iref {
                new_qual = if bmul {
                    (b.qual()[bpos] as f32 * 0.8) as u8
                } else {
                    0
                };

                set_qual(b, bpos, new_qual)?;
                // println! {"Adjusting to deletion in read A: POS: {bpos} QUAL: {new_qual} {}", std::str::from_utf8(a.qname())?}
                if let Some((n_bpos, n_b_iref)) = bp.next() {
                    b_iref = n_b_iref;
                    bpos = n_bpos;
                } else {
                    return Ok(());
                }
            }
        }
        // check for deletion in read B
        // if b_iref > a_iref && bp.passed_deletion() {
        else if b_iref > a_iref && bp.after_del() {
            while a_iref < b_iref {
                new_qual = if amul {
                    (a.qual()[apos] as f32 * 0.8) as u8
                } else {
                    0
                };

                set_qual(a, apos, new_qual)?;
                // println! {"Adjusting to deletion in read B: POS: {apos} QUAL: {new_qual} {}", std::str::from_utf8(a.qname())?}
                if let Some((n_apos, n_a_iref)) = ap.next() {
                    a_iref = n_a_iref;
                    apos = n_apos;
                } else {
                    return Ok(());
                }
            }
        };

        null_ref_bases(
            a,
            apos,
            b,
            bpos,
            amul,
            &mut base_a,
            &mut base_b,
            &mut new_qual,
        )?;
    }

    // Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::record::{Cigar, CigarString};

    #[test]
    pub fn qual_set_test1() {
        let mut record = Record::new();
        record.set(
            b"read1",
            Some(&CigarString(vec![Cigar::Match(4)])),
            b"AAAA",
            b"####",
        );

        set_qual(&mut record, 0, 0).unwrap();
        assert_eq!(record.qual()[0], 0);
    }

    #[test]
    pub fn qual_set_test2() {
        let mut record = Record::new();
        record.set(
            b"read2",
            Some(&CigarString(vec![
                Cigar::Match(4),
                Cigar::Ins(5),
                Cigar::Match(3),
            ])),
            b"AAAAGAAAAAAA",
            b"############",
        );

        set_qual(&mut record, 11, 0).unwrap();
        assert_eq!(record.qual()[11], 0);

        set_qual(&mut record, 8, 4).unwrap();
        assert_eq!(record.qual()[8], 4);

        set_qual(&mut record, 0, 4).unwrap();
        assert_eq!(record.qual()[0], 4);
    }

    #[test]
    pub fn qual_set_bounds_check() {
        let mut record = Record::new();
        record.set(
            b"read2",
            Some(&CigarString(vec![
                Cigar::Match(4),
                Cigar::Ins(5),
                Cigar::Match(3),
            ])),
            b"AAAAGAAAAAAA",
            b"############",
        );

        assert!(set_qual(&mut record, 12, 8).is_err());
    }

    #[test]
    pub fn test_overlap_tweak1() {
        let mut a = Record::new();
        a.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(4),
                Cigar::Del(1),
                Cigar::Match(5),
            ])),
            b"AAAAGTACA",
            b"#########",
        );

        a.set_pos(1);

        let mut b = Record::new();
        b.set(
            b"read0",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"TAAATGTACT",
            b"E########E",
        );

        b.set_pos(1);
        b.set_reverse();

        tweak_overlap_qual(&mut a, &mut b).unwrap();

        let exp_qual_4 = (b'#' as f32 * 0.8) as u8;
        assert_eq!(b.qual()[4], exp_qual_4);
        assert_eq!(a.qual()[4], 0);

        let exp_qual_0 = (b'E' as f32 * 0.8) as u8;
        assert_eq!(b.qual()[0], exp_qual_0);
        assert_eq!(a.qual()[0], 0);

        let exp_qual_8 = (b'E' as f32 * 0.8) as u8;
        assert_eq!(b.qual()[9], exp_qual_8);
        assert_eq!(a.qual()[8], 0);
    }
}

#![allow(dead_code)]
use crate::pileup::PileUp;
use anyhow::Error;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};
use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::slice;
use std::{cell::RefCell, cmp::Ordering, rc::Rc};

pub type OverlapMap = HashMap<u64, Rc<RefCell<PileUp>>>;

pub trait MapOverlaps {
    fn push(&mut self, r: Rc<RefCell<PileUp>>);
    fn delete_hash(&mut self, r: u64);
}

pub fn hash_qname(r: &Record) -> u64 {
    let mut hasher = DefaultHasher::new();
    r.qname().hash(&mut hasher);
    hasher.finish()
}

impl MapOverlaps for OverlapMap {
    fn push(&mut self, plp: Rc<RefCell<PileUp>>) {
        let mut r = &mut plp.borrow_mut().rec;

        if r.is_mate_unmapped() || !r.is_proper_pair() {
            return;
        }

        if r.mtid() >= 0 && (r.mtid() != r.tid()) {
            return;
        }

        let h = hash_qname(&r);

        if let Some(mate) = self.get_mut(&h) {
            tweak_overlap_qual(&mut mate.borrow_mut().rec, &mut r).unwrap();
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
}

/// this is just a modified version of [Record::qual] that returns a mutable slice to the qual
/// array.
pub fn qual_mut(r: &mut Record) -> &mut [u8] {
    let i: usize = r.inner.core.l_qname as usize + r.cigar_len() * 4 + (r.seq_len() + 1) / 2;
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

pub fn tweak_overlap_qual(a: &mut Record, b: &mut Record) -> Result<(), Error> {
    let mut new_qual: u8;
    let mut base_a @ mut base_b: u8;
    let amul @ bmul: bool;

    // println! {"QNAME: {} |", std::str::from_utf8(a.qname())?}

    // we assume that we encounter reads in order (e.g coord-sorted).
    assert!(a.pos() <= b.pos());

    // using aligned_pairs() for this, since I just need sequence that overlaps in covered ref
    // positions.
    // gets iterators over tuples of (read_pos, ref_pos).

    // we use [Record::aligned_pairs_full] to also retrieve deletion bases, for which we also adjust quality
    // (treated the same as mismatches)
    let mut ap = a
        .aligned_pairs_full()
        .map(|x| (x[0].map(|x| x as usize), x[1]))
        .peekable();

    let mut bp = b
        .aligned_pairs_full()
        .map(|x| (x[0].map(|x| x as usize), x[1]))
        .peekable();

    // for now: using htslib's hashing heuristic to decide which read to modify, bound to change
    // (maybe)
    match decide_which_read(a.qname()) {
        true => (amul, bmul) = (true, false),
        false => (amul, bmul) = (false, true),
    }

    loop {
        match (ap.peek(), bp.peek()) {
            (Some((aread, aref)), Some((bread, bref))) => match (aread, aref, bread, bref) {
                (_, None, _, None) => _ = (ap.next(), bp.next()), // both with insertion
                (_, _, _, None) => _ = bp.next(),                 // b has an insertion at this base
                (_, None, _, _) => _ = ap.next(),                 // a has an insertion at this base

                // we are now at a existing reference position for both read A and B
                (aread, Some(aref), bread, Some(bref)) => {
                    if *aref < b.pos() {
                        // we're still behind read B, so advance to catch up, and try again
                        ap.next();
                        continue;
                    }

                    assert_eq!(aref, bref, "{aref} {bref}");

                    // we now are at matching ref positions covered by both reads
                    match (aread, bread) {
                        // read B has a deletion
                        (&Some(aread), None) => {
                            if amul {
                                new_qual = (a.qual()[aread] as f64 * 0.8) as u8;
                            } else {
                                new_qual = 0;
                            }

                            set_qual(a, aread, new_qual)?;
                        }

                        // read A has a deletion
                        (None, &Some(bread)) => {
                            if bmul {
                                new_qual = (b.qual()[bread] as f64 * 0.8) as u8;
                            } else {
                                new_qual = 0;
                            }

                            set_qual(b, bread, new_qual)?;
                        }

                        // both have non-insertion bases at this position
                        (&Some(aread), &Some(bread)) => {
                            (base_a, base_b) = (a.seq()[aread], b.seq()[bread]);

                            // both bases mismatch, so pick based on quality (or random if tie)
                            if base_a != base_b {
                                match a.qual()[aread].cmp(&b.qual()[bread]) {
                                    Ordering::Less => {
                                        new_qual = (b.qual()[bread] as f64 * 0.8) as u8;
                                        set_qual(a, aread, 0)?;
                                        set_qual(b, bread, new_qual)?;
                                    }

                                    Ordering::Greater => {
                                        new_qual = (a.qual()[aread] as f64 * 0.8) as u8;
                                        set_qual(a, aread, new_qual)?;
                                        set_qual(b, bread, 0)?;
                                    }

                                    Ordering::Equal => {
                                        if amul {
                                            new_qual = (a.qual()[aread] as f64 * 0.8) as u8;
                                            set_qual(a, aread, new_qual)?;
                                            set_qual(b, bread, 0)?;
                                        } else {
                                            new_qual = (b.qual()[bread] as f64 * 0.8) as u8;
                                            set_qual(a, aread, 0)?;
                                            set_qual(b, bread, new_qual)?;
                                        }
                                    }
                                }
                            } else {
                                // both bases match; Bump the quality up for one read and null the
                                // other's
                                new_qual = (a.qual()[aread].wrapping_add(b.qual()[bread])).min(200);

                                // set quals accordingly
                                if amul {
                                    set_qual(a, aread, new_qual)?;
                                    set_qual(b, bread, 0)?;
                                } else {
                                    set_qual(a, aread, 0)?;
                                    set_qual(b, bread, new_qual)?;
                                }

                                // println!("Adjusting quality to be ultra-confident {aread} | A POS: {aread} | B POS: {bread} | A QUAL: {} | B QUAL: {}", a.qual()[aread], b.qual()[bread])
                            }
                        }

                        (None, None) => (), // both have dels, so move on
                    }

                    // done with this reference position. Moving on...
                    ap.next();
                    bp.next();
                }
            },

            // once we run out of bases for one of the reads, no use comparing.
            (None, Some(_)) | (Some(_), None) => break,
            (None, None) => break,
        }
    }

    Ok(())
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

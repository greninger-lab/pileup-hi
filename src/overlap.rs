#![allow(dead_code)]
use crate::pileup::PileUp;
use anyhow::Error;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};
use std::collections::HashMap;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::slice;
use std::{cell::RefCell, rc::Rc};

pub type OverlapMap = HashMap<u64, Rc<RefCell<PileUp>>>;

pub trait MapOverlaps {
    fn push(&mut self, r: PileUp) -> OverlapInsResult;
    fn delete_hash(&mut self, r: u64);
}

pub enum OverlapInsResult {
    Inserted(Rc<RefCell<PileUp>>),
    Rejected(PileUp),
}

pub fn hash_qname(r: &Record) -> u64 {
    let mut hasher = DefaultHasher::new();
    r.qname().hash(&mut hasher);
    hasher.finish()
}

impl MapOverlaps for OverlapMap {
    fn push(&mut self, mut plp: PileUp) -> OverlapInsResult {
        let r = &mut plp.rec;
        let h = hash_qname(r);

        if r.is_mate_unmapped() || !r.is_proper_pair() {
            return OverlapInsResult::Rejected(plp);
        }

        if r.mtid() >= 0 && (r.mtid() != r.tid()) {
            return OverlapInsResult::Rejected(plp);
        }

        if let Some(mate) = self.get_mut(&h) {
            tweak_overlap_qual(&mut mate.borrow_mut().rec, r).unwrap();
            return OverlapInsResult::Rejected(plp);
        }

        if r.pos() < r.mpos() || (r.is_paired() && r.mpos() == -1) {
            // criteria passed, insert
            plp.in_overlap = true;
            let ins = Rc::new(RefCell::new(plp));
            let _ = self.insert(h, Rc::clone(&ins));
            OverlapInsResult::Inserted(ins)
        } else {
            OverlapInsResult::Rejected(plp)
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
    let mut chars = chars.into_iter();
    let mut h = match chars.next() {
        Some(c) => *c as u64,
        None => 0,
    };
    for c in chars {
        h = (h << 5).wrapping_sub(h) + *c as u64;
    }
    h = h.wrapping_add(!(h << 15));
    h ^= h >> 10;
    h = h.wrapping_add(h << 3);
    h ^= h >> 6;
    h = h.wrapping_add(!(h << 11));
    h ^= h >> 16;
    h & 1 != 0
}

pub fn tweak_overlap_qual(a: &mut Record, b: &mut Record) -> Result<(), Error> {
    let mut new_qual: u8;
    let mut base_a @ mut base_b: u8;
    let amul @ bmul: bool;

    // using aligned_pairs() for this, since I just need sequence that overlaps in covered ref
    // positions.
    // gets iterators over tuples of (read_pos, ref_pos).

    // we use [Record::aligned_pairs_full] to also retrieve deletion bases, for which we also adjust quality
    // (treated the same as mismatches)
    let ap = a
        .aligned_pairs_full()
        .map(|x| (x[0].map(|x| x as usize), x[1]));

    let bp = b.aligned_pairs();

    // for now: using htslib's hashing heuristic to decide which read to modify, bound to change
    // (maybe)
    match decide_which_read(a.qname()) {
        true => (amul, bmul) = (true, false),
        false => (amul, bmul) = (false, true),
    }

    let bhash: HashMap<i64, usize> =
        HashMap::from_iter(bp.into_iter().map(|x| (x[1], x[0] as usize)));

    for (read_pos, ref_pos) in ap {
        if let Some(rp) = ref_pos {
            // read a has a deletion relative to read b
            if read_pos.is_none() && bhash.contains_key(&rp) {
                let bpos = *bhash.get(&rp).unwrap();

                // set qual accordingly
                if bmul {
                    new_qual = (b.qual()[bpos] as f32 * 0.8) as u8;
                } else {
                    new_qual = 0;
                }

                set_qual(b, bpos, new_qual)?;
                continue;
            }
            // read A has a base at this ref position
            if let Some(read_pos) = read_pos {
                // both reads have bases at the given reference position
                if let Some(other_read_pos) = bhash.get(&rp) {
                    let (read_pos, other_read_pos) = (read_pos, *other_read_pos);
                    (base_a, base_b) = (a.seq()[read_pos], b.seq()[other_read_pos]);

                    // both bases match :)
                    if base_a == base_b {
                        new_qual = a.qual()[read_pos].wrapping_add(b.qual()[other_read_pos]);

                        // set quals accordingly
                        if amul {
                            set_qual(a, read_pos, new_qual)?;
                            set_qual(b, other_read_pos, 0)?;
                        } else {
                            set_qual(a, read_pos, 0)?;
                            set_qual(b, other_read_pos, new_qual)?;
                        }
                    } else {
                        // bases are a mismatch :(
                        match a.qual()[read_pos].cmp(&b.qual()[other_read_pos]) {
                            // higher confidence in A
                            std::cmp::Ordering::Greater => {
                                new_qual = (a.qual()[read_pos] as f32 * 0.8) as u8;
                                set_qual(a, read_pos, new_qual)?;
                                set_qual(b, other_read_pos, 0)?;
                            }

                            // higher confidence in B
                            std::cmp::Ordering::Less => {
                                new_qual = (b.qual()[other_read_pos] as f32 * 0.8) as u8;
                                set_qual(a, read_pos, 0)?;
                                set_qual(b, other_read_pos, new_qual)?;
                            }

                            // equal score, so pick based on hash
                            std::cmp::Ordering::Equal => {
                                if amul {
                                    new_qual = (a.qual()[read_pos] as f32 * 0.8) as u8;
                                    set_qual(a, read_pos, new_qual)?;
                                    set_qual(b, other_read_pos, 0)?;
                                } else {
                                    new_qual = (b.qual()[other_read_pos] as f32 * 0.8) as u8;
                                    set_qual(a, read_pos, 0)?;
                                    set_qual(b, other_read_pos, new_qual)?;
                                }
                            }
                        }
                    }
                } else {
                    // read b has a deletion relative to read a
                    if let Some(_) = bhash.get(&ref_pos.unwrap()) {
                        let rp = read_pos;
                        if amul {
                            new_qual = (a.qual()[rp] as f32 * 0.8) as u8;
                        } else {
                            new_qual = 0;
                        }

                        set_qual(b, read_pos, new_qual)?;
                        continue;
                    }
                }
            }
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
}

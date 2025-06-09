#![allow(dead_code)]
use crate::pileup::PileUp;
use anyhow::Error;
use rust_htslib::bam::{record::Cigar, Record};
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

/// a twist on [rust_htslib::bam::ext::IterAlignedPairsFull] that differentiates between
/// [Cigar::RefSkip] and [Cigar::Del].
pub struct CigarWalker {
    del_remaining: u32,
    ins_remaining: u32,
    match_remaining: u32,
    refskip_remaining: u32,
    ref_pos: i64,
    read_pos: usize,
    cigar: Vec<Cigar>,
    cigar_index: usize,
    in_del: bool,
}

impl CigarWalker {
    pub fn new(r: &Record) -> Self {
        Self {
            del_remaining: 0,
            ins_remaining: 0,
            match_remaining: 0,
            refskip_remaining: 0,
            ref_pos: r.pos() - 1,
            read_pos: usize::MAX,
            cigar: r.cigar().take().0,
            cigar_index: 0,
            in_del: false,
        }
    }

    pub fn move_to_next_match(&mut self) -> Option<Cigar> {
        // if self.match_remaining > 0 {
        //     self.match_remaining -= 1;
        //     self.read_pos -= 1;
        //     self.

        //     return Some(Cigar::Match(1)); // already in match block
        // }

        let mut ret: Option<Cigar>;
        loop {
            ret = self.next();
            match ret {
                None => {
                    return None;
                }

                Some(Cigar::Match(_)) => {
                    return ret;
                }

                _ => continue,
            }
        }
    }

    pub fn move_to_ref_pos(&mut self, pos: i64) -> Option<Cigar> {
        if pos == self.ref_pos {
            return Some(Cigar::Match(1));
        }

        let mut ret = self.next();
        while self.ref_pos < pos {
            ret = self.next();
            match ret {
                None => break,
                _ => continue,
            }
        }
        ret
    }
}

impl Iterator for CigarWalker {
    type Item = Cigar;
    fn next(&mut self) -> Option<Self::Item> {
        if self.del_remaining > 0 {
            self.del_remaining -= 1;
            self.ref_pos += 1;
            return Some(Cigar::Del(1));
        }
        if self.ins_remaining > 0 {
            self.ins_remaining -= 1;
            self.read_pos += 1;
            return Some(Cigar::Ins(1));
        }
        if self.match_remaining > 0 {
            self.match_remaining -= 1;
            self.ref_pos += 1;
            self.read_pos += 1;
            return Some(Cigar::Match(1));
        }
        if self.refskip_remaining > 0 {
            self.refskip_remaining -= 1;
            self.ref_pos += 1;
            return Some(Cigar::RefSkip(1));
        }

        while self.cigar_index < self.cigar.len() {
            let entry = self.cigar[self.cigar_index];
            match entry {
                Cigar::Match(len) | Cigar::Equal(len) | Cigar::Diff(len) => {
                    self.in_del = false;
                    // self.read_pos += 1;
                    self.read_pos = self.read_pos.wrapping_add(1);
                    self.ref_pos += 1;
                    self.match_remaining = len - 1;
                    self.cigar_index += 1;
                    return Some(Cigar::Match(1));
                }

                Cigar::Ins(len) | Cigar::SoftClip(len) => {
                    self.in_del = false;
                    self.read_pos += 1;
                    self.ins_remaining = len - 1;
                    self.cigar_index += 1;
                    return Some(Cigar::Ins(1));
                }

                Cigar::Del(len) => {
                    self.in_del = true;
                    self.ref_pos += 1;
                    self.del_remaining = len - 1;
                    self.cigar_index += 1;
                    return Some(Cigar::Del(1));
                }

                Cigar::RefSkip(len) => {
                    self.in_del = false;
                    self.ref_pos += 1;
                    self.refskip_remaining = len - 1;
                    self.cigar_index += 1;
                    return Some(Cigar::RefSkip(1));
                }

                Cigar::HardClip(_) => (),
                Cigar::Pad(_) => panic!("Unsupported padding op"),
            }
            self.cigar_index += 1
        }
        None
    }
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
            let _ = self.insert(h, Rc::clone(&plp));
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
    let chars = chars.into_iter();
    let mut h = 0u32;
    for c in chars {
        h = (h << 5) - (h) + *c as u32;
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
    let mut ret_a @ mut ret_b: Option<Cigar>;
    // println! {"QNAME: {} |", std::str::from_utf8(a.qname())?}

    // we assume that we encounter reads in order (e.g coord-sorted).
    assert!(a.pos() <= b.pos());
    // println! {"==============================="}

    let mut a_iter = CigarWalker::new(&a);
    let mut b_iter = CigarWalker::new(&b);

    let (amul, bmul) = match decide_which_read(a.qname()) {
        true => (true, false),
        false => (false, true),
    };

    while a_iter.ref_pos != b_iter.ref_pos {
        match a_iter.next() {
            None => break,
            _ => (),
        };
    }

    loop {
        ret_a = a_iter.move_to_next_match();
        // println! {"A: {}", a_iter.read_pos}
        match ret_a {
            None => break,
            _ => (),
        }

        ret_b = b_iter.move_to_next_match();
        // println! {"B: {}", b_iter.read_pos}
        match ret_b {
            None => break,
            _ => (),
        }

        // println! {"{} {}", a_iter.ref_pos, b_iter.ref_pos}

        if a_iter.ref_pos != b_iter.ref_pos {
            // del in read B
            if a_iter.ref_pos < b_iter.ref_pos {
                while a_iter.ref_pos < b_iter.ref_pos {
                    let new_qual = match amul {
                        true => (a.qual()[a_iter.read_pos] as f32 * 0.8) as u8,
                        false => 0,
                    };
                    // println! {"del in read b: B: {}", b_iter.read_pos}
                    set_qual(a, a_iter.read_pos, new_qual)?;
                    // println! {"adjusted qual at A {} to {}", a_iter.read_pos, new_qual}
                    if !matches!(a_iter.next(), Some(Cigar::Match(_))) {
                        break;
                    }
                    // if a_iter.next().is_none() {
                    //     break;
                    // }
                }
                // a_iter.next();
            }

            // del in read A
            if b_iter.ref_pos < a_iter.ref_pos {
                while b_iter.ref_pos < a_iter.ref_pos {
                    new_qual = match bmul {
                        true => (b.qual()[b_iter.read_pos] as f32 * 0.8) as u8,
                        false => 0,
                    };
                    // println! {"del in read a: A: {}", a_iter.read_pos}
                    set_qual(b, b_iter.read_pos, new_qual)?;
                    // println! {"adjusted qual at B {} to {}", b_iter.read_pos, new_qual}
                    if !matches!(b_iter.next(), Some(Cigar::Match(_))) {
                        break;
                    }
                }
                // b_iter.next();
            }
        }

        if a_iter.ref_pos != b_iter.ref_pos {
            continue;
        }

        // println! {"POS AFTER DELETION CORRECTION: {} {}", a_iter.ref_pos, b_iter.ref_pos}
        assert_eq!(
            a_iter.ref_pos,
            b_iter.ref_pos,
            "{}->{} {}->{}",
            a.pos(),
            a.cigar(),
            b.pos(),
            b.cigar()
        );

        // read does not match
        if a.seq()[a_iter.read_pos] != b.seq()[b_iter.read_pos] {
            // println! {"mismatch REF: {} {} READ: {} {}", a_iter.ref_pos, b_iter.ref_pos, a_iter.read_pos, b_iter.read_pos}
            match a.qual()[a_iter.read_pos].cmp(&b.qual()[b_iter.read_pos]) {
                Ordering::Greater => {
                    new_qual = (a.qual()[a_iter.read_pos] as f32 * 0.8) as u8;
                    set_qual(a, a_iter.read_pos, new_qual)?;
                    set_qual(b, b_iter.read_pos, 0)?;
                }
                Ordering::Less => {
                    new_qual = (b.qual()[b_iter.read_pos] as f32 * 0.8) as u8;
                    set_qual(a, a_iter.read_pos, 0)?;
                    set_qual(b, b_iter.read_pos, new_qual)?;
                }
                Ordering::Equal => match amul {
                    true => {
                        new_qual = (a.qual()[a_iter.read_pos] as f32 * 0.8) as u8;
                        set_qual(a, a_iter.read_pos, new_qual)?;
                        set_qual(b, b_iter.read_pos, 0)?;
                    }

                    false => {
                        new_qual = (b.qual()[b_iter.read_pos] as f32 * 0.8) as u8;
                        set_qual(a, a_iter.read_pos, 0)?;
                        set_qual(b, b_iter.read_pos, new_qual)?;
                    }
                },
            }
        } else {
            // println! {"match REF: {} {} READ: {} {}", a_iter.ref_pos, b_iter.ref_pos, a_iter.read_pos, b_iter.read_pos}
            new_qual = (a.qual()[a_iter.read_pos] + b.qual()[b_iter.read_pos]).min(200);
            match amul {
                true => {
                    set_qual(a, a_iter.read_pos, new_qual)?;
                    set_qual(b, b_iter.read_pos, 0)?;
                }

                false => {
                    set_qual(a, a_iter.read_pos, 0)?;
                    set_qual(b, b_iter.read_pos, new_qual)?;
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

    const E_ADJ: u8 = (b'E' as f32 * 0.8) as u8;
    const E_CONF: u8 = b'E' * 2;
    const HASH_ADJ: u8 = (b'#' as f32 * 0.8) as u8;

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
            b"read1",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"TAAATGTACT",
            b"E########E",
        );

        b.set_pos(1);
        b.set_reverse();

        tweak_overlap_qual(&mut a, &mut b).unwrap();

        // let exp_qual_4 = (b'#' as f32 * 0.8) as u8;
        assert_eq!(a.qual()[4], 0, "{} {}", a.qual()[4], b.qual()[4]);
        assert_eq!(b.qual()[4], HASH_ADJ);

        // let exp_qual_0 = (b'E' as f32 * 0.8) as u8;
        assert_eq!(a.qual()[0], 0);
        assert_eq!(b.qual()[0], E_ADJ);

        // let exp_qual_8 = (b'E' as f32 * 0.8) as u8;
        assert_eq!(b.qual()[9], E_ADJ);
        assert_eq!(a.qual()[8], 0);
    }

    #[test]
    pub fn test_overlap_tweak2() {
        let mut a = Record::new();
        a.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(4),
                Cigar::Del(2),
                Cigar::Match(5),
            ])),
            b"AAAATACA",
            b"########",
        );

        a.set_pos(1);

        let mut b = Record::new();

        b.set(
            b"read1",
            Some(&CigarString(vec![Cigar::Match(10)])),
            b"TAAATGTACT",
            b"E########E",
        );

        b.set_pos(1);
        b.set_reverse();

        tweak_overlap_qual(&mut a, &mut b).unwrap();

        // check deletion here...
        // let exp_qual_4 = (b'#' as f32 * 0.8) as u8;
        assert_eq!(a.qual()[4], 0);
        assert_eq!(b.qual()[4], HASH_ADJ);

        // let exp_qual_0 = (b'E' as f32 * 0.8) as u8;
        assert_eq!(a.qual()[0], 0);
        assert_eq!(b.qual()[0], E_ADJ);

        // let exp_qual_8 = (b'E' as f32 * 0.8) as u8;
        assert_eq!(b.qual()[9], E_ADJ);
        assert_eq!(a.qual()[7], 0);
    }

    #[test]
    // now test with ref offset of 2
    pub fn test_overlap_tweak3() {
        let mut a = Record::new();

        a.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(6),
                Cigar::Del(1),
                Cigar::Match(4),
            ])),
            b"GCTGCAGTAT",
            b"EEEEEEEEEE",
            // b"read1",
            // Some(&CigarString(vec![
            //     Cigar::Match(2),
            //     Cigar::Del(1),
            //     Cigar::Match(3),
            //     Cigar::Del(2),
            //     Cigar::Match(1),
            // ])),
            // b"TGAGGT",
            // b"EEEEEE",
        );

        a.set_pos(1);

        let mut b = Record::new();
        b.set(
            b"read1",
            Some(&CigarString(vec![
                Cigar::Match(2),
                Cigar::Del(1),
                Cigar::Match(3),
                Cigar::Del(2),
                Cigar::Match(1),
            ])),
            b"TGAGGT",
            b"EEEEEE",
        );

        b.set_pos(3);

        tweak_overlap_qual(&mut a, &mut b).unwrap();

        // check deletion here...
        assert_eq!(b.qual()[0], E_CONF);
        assert_eq!(a.qual()[0], b'E');

        // first set of deletions
        assert_eq!(a.qual()[3], E_ADJ);
        assert_eq!(b.qual()[4], E_CONF);

        // second set of deletions
        assert_eq!(b.qual()[8], E_ADJ);
        assert_eq!(b.qual()[9], E_ADJ);
    }
}

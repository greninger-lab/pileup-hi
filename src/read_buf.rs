use crate::alignment::{cigar2rlen, CigarState, PileupAlignment, PileupAlignmentRef, CIGAR_STATE_UNINIT};
use crate::overlap::{MapOverlaps, OverlapMap};
use rust_htslib::bam::Record;
use std::{cell::RefCell, collections::HashMap, rc::Rc};

pub struct ReadBuffer {
    pub rbuf: Vec<PileupAlignmentRef>,
    pub len: i64,
    pub backup_buf: Vec<PileupAlignmentRef>,
    pub overlap_map: Option<OverlapMap>,
    pub depth: usize,
    pub max_depth: usize,
}

pub enum BufPushResult {
    Pushed(PileupAlignmentRef),
    DifferentReference(PileupAlignmentRef),
    Unmapped,
    MaxDepthMet,
    BeforePos,
}

impl ReadBuffer {
    #[inline(always)]
    pub fn attempt_push(&mut self, r: &Record, pos: i64, tid: i32) -> BufPushResult {
        let mut dif_ref = false;

        if r.is_unmapped() {
            return BufPushResult::Unmapped;
        }

        if r.tid() != tid {
            dif_ref = true;
        }

        let read_len_from_cigar = cigar2rlen(r);

        if read_len_from_cigar > self.len {
            self.len = read_len_from_cigar;
        }

        if !dif_ref && r.pos() + read_len_from_cigar - 1 < pos {
            return BufPushResult::BeforePos;
        }

        if r.tid() == tid && r.pos() == pos && self.depth >= self.max_depth {
            if let Some(ov) = &mut self.overlap_map {
                ov.delete_read(r);
            }
            return BufPushResult::MaxDepthMet;
        }

        let cstate = CigarState {
            cig: r.cigar(),
            icig: CIGAR_STATE_UNINIT,
            iseq: 0,
            bam_pos: r.pos(),
            read_len_from_cigar,
        };

        let plp = PileupAlignment::new(r.clone(), cstate);

        let plp_ref = Rc::new(RefCell::new(plp));

        if let Some(overlap_map) = &mut self.overlap_map {
            overlap_map.push(Rc::clone(&plp_ref));
        }

        self.rbuf.push(Rc::clone(&plp_ref));
        self.depth += 1;

        if dif_ref {
            BufPushResult::DifferentReference(plp_ref)
        } else {
            BufPushResult::Pushed(plp_ref)
        }
    }

    pub fn new(depth: usize, disable_overlaps: bool) -> Self {
        let rbuf: Vec<PileupAlignmentRef> = Vec::with_capacity(500);
        let backup_buf: Vec<PileupAlignmentRef> = Vec::with_capacity(500);

        let max_depth = if depth.cmp(&0).is_eq() { usize::MAX } else { depth };
        let len = 0;

        let overlap_map = match disable_overlaps {
            false => Some(HashMap::new()),
            true => None,
        };

        Self {
            rbuf,
            backup_buf,
            overlap_map,
            len,
            depth: 0,
            max_depth,
        }
    }

    pub fn start(&self) -> i64 {
        if let Some(r) = self.rbuf.first() {
            r.borrow().rec.pos()
        } else {
            i64::MAX
        }
    }

    pub fn reset(&mut self) {
        assert!(self.rbuf.is_empty());
        std::mem::swap(&mut self.rbuf, &mut self.backup_buf);
        // if let Some(ov) = &mut self.overlap_map {
        //     ov.shrink_to_fit();
        // }
    }
}

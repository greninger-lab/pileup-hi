#![allow(dead_code)]
use crate::overlap::{MapOverlaps, OverlapMap};
use crate::pileup::{CigarState, PileUp, cigar2rlen};
use rust_htslib::bam::Record;
use std::{cell::RefCell, collections::HashMap, rc::Rc};

pub struct ReadBuffer {
    pub rbuf: Vec<Rc<RefCell<PileUp>>>,
    pub len: i64,
    pub backup_buf: Vec<Rc<RefCell<PileUp>>>,
    pub overlap_map: Option<OverlapMap>,
    pub depth: usize,
    pub max_depth: usize,
}

#[derive(Debug, Eq, PartialEq)]
pub enum BufPushResult {
    // AfterWindow(usize),
    Pushed,
    DifferentReference,
    Unmapped,
    MaxDepthMet,
}

impl ReadBuffer {
    pub fn attempt_push(&mut self, r: &Record, pos: i64, tid: i32) -> BufPushResult {
        if r.is_unmapped() {
            return BufPushResult::Unmapped;
        }

        if r.tid() != tid {
            return BufPushResult::DifferentReference;
        }

        if cigar2rlen(&r) > self.len {
            self.len = cigar2rlen(&r);
        }

        if r.pos() + self.len - 1 < pos {
            panic!(); // unsorted
        }

        if r.tid() == tid && r.pos() == pos && self.depth >= self.max_depth {
            if let Some(ov) = &mut self.overlap_map {
                ov.delete_read(r);
            }
            return BufPushResult::MaxDepthMet;
        }

        let cstate = CigarState {
            cig: r.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: r.pos() as u32,
            qpos: 0,
            del: false,
        };

        let plp = PileUp {
            rec: r.clone(),
            cstate,
        };

        let plp_ref = Rc::new(RefCell::new(plp));

        if let Some(overlap_map) = &mut self.overlap_map {
            overlap_map.push(Rc::clone(&plp_ref));
        }

        self.rbuf.push(plp_ref);
        self.depth += 1;

        BufPushResult::Pushed
    }

    pub fn new(depth: usize, disable_overlaps: bool) -> Self {
        let rbuf: Vec<Rc<RefCell<PileUp>>> = Vec::with_capacity(500);
        let backup_buf: Vec<Rc<RefCell<PileUp>>> = Vec::with_capacity(500);
        let max_depth = depth.cmp(&0).is_eq().then_some(usize::MAX).unwrap_or(depth);
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

    pub fn reset(&mut self) {
        assert!(self.rbuf.is_empty());
        std::mem::swap(&mut self.rbuf, &mut self.backup_buf);
    }
}

#![allow(dead_code)]
use crate::overlap::{MapOverlaps, OverlapMap};
use crate::pileup::{cigar2rlen, CigarState, PileUp};
use rust_htslib::bam::Record;
use std::{cell::RefCell, collections::HashMap, rc::Rc};

pub struct ReadBuffer {
    pub rbuf: Vec<Rc<RefCell<PileUp>>>,
    pub len: usize,
    pub backup_buf: Vec<Rc<RefCell<PileUp>>>,
    pub overlap_map: OverlapMap,
}

#[derive(Debug, Eq, PartialEq)]
pub enum BufPushResult {
    AfterWindow(usize),
    Pushed,
    DifferentReference,
    Unmapped,
}

impl ReadBuffer {
    pub fn c_to_next_window(&mut self, next_pos: i64, cur_pos: usize) -> usize {
        let next_pos = next_pos as usize;
        std::cmp::max(0, next_pos - (cur_pos + self.len - 1))
    }

    pub fn attempt_push(&mut self, r: &Record, pos: usize, tid: u32) -> BufPushResult {
        if r.is_unmapped() {
            return BufPushResult::Unmapped;
        }

        if r.tid() as u32 != tid {
            return BufPushResult::DifferentReference;
        }

        if cigar2rlen(&r) > self.len {
            self.len = cigar2rlen(&r);
        }

        if r.pos() as usize + self.len - 1 < pos {
            panic!(); // unsorted
        }

        if r.pos() as usize > pos + self.len - 1 {
            let window_start = self.c_to_next_window(r.pos(), pos);
            return BufPushResult::AfterWindow(window_start);
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

        self.overlap_map.push(Rc::clone(&plp_ref));
        self.rbuf.push(plp_ref);

        BufPushResult::Pushed
    }

    pub fn new() -> Self {
        let rbuf: Vec<Rc<RefCell<PileUp>>> = Vec::with_capacity(500);
        let backup_buf: Vec<Rc<RefCell<PileUp>>> = Vec::with_capacity(500);
        let overlap_map = HashMap::new();
        let len = 0;

        Self {
            rbuf,
            backup_buf,
            overlap_map,
            len,
        }
    }

    pub fn reset(&mut self) {
        assert!(self.rbuf.is_empty());
        std::mem::swap(&mut self.rbuf, &mut self.backup_buf);
    }
}

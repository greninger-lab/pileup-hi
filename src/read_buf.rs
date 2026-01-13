use crate::alignment::{cigar2rlen, CigarState, PileupAlignment, PileupAlignmentRef, CIGAR_STATE_UNINIT};
use crate::cigar_resolve::resolve_cigar;
use crate::overlap::{MapOverlaps, OverlapMap};
use anyhow::Error;
use log::error;
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
    Pushed,
    DifferentReference,
    Unmapped,
    MaxDepthMet,
    BeforePos,
}

impl ReadBuffer {
    #[inline(always)]
    pub fn store(&mut self, r: &Record, read_len_from_cigar: i64, pos: i64) {
        let cstate = CigarState {
            cig: r.cigar(),
            icig: CIGAR_STATE_UNINIT,
            iseq: 0,
            bam_pos: r.pos(),
            read_len_from_cigar,
        };

        let mut plp = PileupAlignment::new(r.clone(), cstate);

        // for the love of god, don't remove this.
        //
        // if we encounter a read that starts before pos, we need to initialize its cigar state to
        // catch up to the current position, as resolve_cigar increments in one-cigar-op intervals.
        //
        // if we don't, and if the query position is at cigar block 2+, the cigar state will stay
        // stale.
        if plp.rec.pos() < pos {
            for i in plp.rec.pos()..pos {
                resolve_cigar(&mut plp, i);
            }
        }

        let plp_ref = Rc::new(RefCell::new(plp));

        if let Some(overlap_map) = &mut self.overlap_map {
            overlap_map.push(Rc::clone(&plp_ref));
        }

        self.rbuf.push(Rc::clone(&plp_ref));
        self.depth += 1;
    }

    #[inline(always)]
    pub fn attempt_push(&mut self, tid: i32, pos: i64, r: &Record) -> Result<BufPushResult, Error> {
        if r.is_unmapped() {
            return Ok(BufPushResult::Unmapped);
        }

        // check for unsorted input
        if let Some((head_tid, head_pos)) = self.tail() {
            if r.tid() < head_tid {
                error!(
                    "File unsorted by reference: tid {} comes after tid {}. Read name: {}",
                    r.tid(),
                    head_tid,
                    std::str::from_utf8(r.qname()).unwrap()
                );

                anyhow::bail!("Unsorted");
            }

            if r.pos() < head_pos && head_tid == r.tid() {
                error!(
                    "File unsorted by coordinate: pos {} comes after pos {}. Read name: {}",
                    r.pos(),
                    head_pos,
                    std::str::from_utf8(r.qname()).unwrap()
                );
                anyhow::bail!("Unsorted");
            }
        }

        let read_len_from_cigar = cigar2rlen(r);

        if read_len_from_cigar > self.len {
            self.len = read_len_from_cigar;
        }

        if r.tid() != tid {
            self.store(r, read_len_from_cigar, pos);
            return Ok(BufPushResult::DifferentReference);
        }

        if r.pos() == pos && self.depth >= self.max_depth {
            if let Some(ov) = &mut self.overlap_map {
                ov.delete_read(r);
            }
            return Ok(BufPushResult::MaxDepthMet);
        }

        if r.pos() + read_len_from_cigar - 1 < pos {
            return Ok(BufPushResult::BeforePos);
        }

        self.store(r, read_len_from_cigar, pos);
        Ok(BufPushResult::Pushed)
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

    pub fn head(&self) -> Option<(i32, i64)> {
        self.rbuf.first().map(|p| {
            let r = &p.borrow().rec;
            (r.tid(), r.pos())
        })
    }

    pub fn tail(&self) -> Option<(i32, i64)> {
        self.rbuf.last().map(|p| {
            let r = &p.borrow().rec;
            (r.tid(), r.pos())
        })
    }

    pub fn reset(&mut self) {
        assert!(self.rbuf.is_empty());
        std::mem::swap(&mut self.rbuf, &mut self.backup_buf);
    }
}

use rust_htslib::bam::{record::CigarStringView, Record};

pub struct PileUpRecord {
    pub rec: Record,
    pub cstate: CigarState,
}

pub struct ReadBuffer {
    pub rbuf: Vec<PileUpRecord>,
    pub len: usize,
}

#[derive(Debug, Eq, PartialEq)]
pub enum BufPushResult {
    BeforeWindow,
    AfterWindow,
    Pushed,
    DifferentReference,
}

pub struct CigarState {
    pub cig: CigarStringView,
    pub icig: usize,  // position in cigar string
    pub iseq: u32,    // position in read sequence that corresponds to cigar pos
    pub bam_pos: u32, // ref coord of first base
}

impl ReadBuffer {
    pub fn push(&mut self, r: Record, pos: usize, tid: u32) -> BufPushResult {
        if r.tid() as u32 != tid {
            return BufPushResult::DifferentReference;
        }

        if r.seq_len() > self.len {
            self.len = r.seq_len();
        }

        if r.pos() as usize + self.len < pos {
            return BufPushResult::BeforeWindow;
        }

        if r.pos() as usize > pos + self.len {
            return BufPushResult::AfterWindow;
        }

        let cstate = CigarState {
            cig: r.cigar(),
            icig: 0,
            iseq: 0,
            bam_pos: r.pos() as u32,
        };

        self.rbuf.push(PileUpRecord { rec: r, cstate });
        BufPushResult::Pushed
    }

    pub fn new() -> Self {
        let rbuf: Vec<PileUpRecord> = Vec::with_capacity(500);
        let len = 0;

        Self { rbuf, len }
    }
}

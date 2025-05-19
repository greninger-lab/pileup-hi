use rust_htslib::bam::Record;

pub struct ReadBuffer {
    pub rbuf: Vec<Record>,
    pub len: usize,
    pub pos: usize,
    pub tid: u32,
}

pub enum BufPushResult {
    BeforeWindow,
    AfterWindow,
    Pushed,
    DifferentReference,
}

impl ReadBuffer {
    pub fn push(&mut self, r: Record) -> BufPushResult {
        if r.tid() as u32 != self.tid {
            return BufPushResult::DifferentReference;
        }

        if r.seq_len() > self.len {
            self.len = r.seq_len();
        }

        if r.pos() as usize + self.len < self.pos {
            return BufPushResult::BeforeWindow;
        }

        if r.pos() as usize > self.pos + self.len {
            return BufPushResult::AfterWindow;
        }

        self.rbuf.push(r);
        BufPushResult::Pushed
    }

    pub fn new() -> Self {
        let rbuf: Vec<Record> = Vec::with_capacity(500);
        let len = 0;
        let pos = 0;
        let tid = 0;

        Self {
            rbuf,
            len,
            pos,
            tid,
        }
    }

    /// Process all reads in a rbuf (analogous to mplp_set_pileup)
    /// https://github.com/samtools/bcftools/blob/05621cfee236a5826d68263d6f566be1443be717/mpileup2/mpileup.c#L880
    pub fn gen_pileup(&mut self) {
        // remember to remove reads here that are at the very end of their c-state
        //
        for r in self.rbuf.iter() {
            println! {"POS: {}, LEN: {}, CIG: {:?}", r.pos(), r.seq_len(), r.cigar()}
        }
        unimplemented!();
    }
}

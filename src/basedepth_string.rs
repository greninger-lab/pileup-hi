use crate::{alignment::PileupAlignment, output::OrderedPileupOutput};
use anyhow::Error;
use indexmap::IndexMap;
use rust_htslib::bam::record::Cigar;
use std::{io::Write, ops::AddAssign};

pub struct BaseDepthString {
    tid: i32,
    pos: i64,
    ref_name: String,
    depth: u32,
    a: u32,
    g: u32,
    c: u32,
    t: u32,
    n: u32,
    gap: u32,
    refskip: u32,
    insertions: IndexMap<Vec<u8>, u32>,
    deletions: IndexMap<Vec<u8>, u32>,
    lock: std::io::StdoutLock<'static>,
}

impl OrderedPileupOutput for BaseDepthString {
    fn tid(&self) -> i32 {
        self.tid
    }

    fn pos(&self) -> i64 {
        self.pos
    }

    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, _ref_seq: Option<&[u8]>) {
        self.update(tid, pos, ref_name);
    }

    fn intake(&mut self, p: &PileupAlignment, refseq: Option<&[u8]>) -> Result<(), Error> {
        self.intake(p, refseq)
    }

    fn write(&mut self) -> Result<(), Error> {
        self.write()
    }

    fn depth(&self) -> u32 {
        self.depth
    }
}

impl BaseDepthString {
    pub fn new() -> Self {
        let s = std::io::stdout();
        let lock = s.lock();
        Self {
            lock,
            tid: 0,
            pos: 0,
            ref_name: "".to_string(),
            a: 0,
            g: 0,
            c: 0,
            t: 0,
            n: 0,
            depth: 0,
            gap: 0,
            refskip: 0,
            insertions: IndexMap::new(),
            deletions: IndexMap::new(),
        }
    }

    fn reset(&mut self) {
        self.a = 0;
        self.g = 0;
        self.c = 0;
        self.t = 0;
        self.n = 0;
        self.depth = 0;
        self.gap = 0;
        self.refskip = 0;
        self.insertions.clear();
        self.deletions.clear();
    }

    pub fn write(&mut self) -> Result<(), Error> {
        let mut buf = itoa::Buffer::new();

        self.lock.write_all(self.ref_name.as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.pos + 1).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.depth).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.a).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.g).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.c).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.t).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.n).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.gap).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(buf.format(self.refskip).as_bytes())?;
        self.lock.write_all(b"\t")?;

        self.lock.write_all(b"[")?;

        let n = self.insertions.len() - 1;
        for (i, (ins, cnt)) in self.insertions.iter().enumerate() {
            self.lock.write_all(buf.format(*cnt).as_bytes())?;
            self.lock.write_all(ins)?;
            if i < n {
                self.lock.write_all(b" ")?
            };
        }

        self.lock.write_all(b"]\t[")?;

        let n = self.deletions.len() - 1;
        for (i, (del, cnt)) in self.deletions.iter().enumerate() {
            self.lock.write_all(buf.format(*cnt).as_bytes())?;
            self.lock.write_all(del)?;
            if i < n {
                self.lock.write_all(b" ")?
            };
        }

        write!(self.lock, "]")?;
        writeln!(self.lock)?;

        self.reset();

        Ok(())
    }

    pub fn update(&mut self, tid: i32, ref_pos: i64, ref_name: &str) {
        self.tid = tid;
        self.pos = ref_pos;
        if self.ref_name != ref_name {
            self.ref_name = ref_name.to_string()
        }
    }

    pub fn intake(&mut self, p: &PileupAlignment, refseq: Option<&[u8]>) -> Result<(), Error> {
        self.depth += 1;
        self.register_pileup(p, refseq)
    }

    pub fn register_pileup(&mut self, p: &PileupAlignment, refseq: Option<&[u8]>) -> Result<(), Error> {
        match p.del {
            false => {
                let readbase = p.rec.seq()[p.qpos];

                match readbase.to_ascii_uppercase() {
                    b'A' => self.a += 1,
                    b'G' => self.g += 1,
                    b'C' => self.c += 1,
                    b'T' => self.t += 1,
                    b'N' => self.n += 1,
                    other => anyhow::bail!("Unrecognized nucleotide character {}", other),
                }
            }

            true => {
                if p.refskip {
                    self.refskip += 1;
                } else {
                    self.gap += 1;
                }
            }
        }

        let mut del_len = -p.indel;

        if p.indel > 0 {
            let mut temp_buf: Vec<u8> = Vec::with_capacity(p.indel as usize);
            expand_insertions(p, &mut temp_buf, &mut del_len)?;
            self.insertions.entry(temp_buf).or_insert(0).add_assign(1);
        }

        if del_len > 0 {
            let mut temp_buf: Vec<u8> = Vec::with_capacity(del_len as usize);
            let mut refbase;

            for i in 1..=del_len as usize {
                refbase = if let Some(refseq) = refseq {
                    refseq[self.pos as usize + i]
                } else {
                    b'N'
                };
                temp_buf.push(refbase);
            }

            self.deletions.entry(temp_buf).or_insert(0).add_assign(1);
        }
        Ok(())
    }
}

pub fn expand_insertions(p: &PileupAlignment, seq_buf: &mut Vec<u8>, ndel: &mut i32) -> Result<(), Error> {
    let mut read_pos: usize;
    let mut read_base: u8;
    let ncig = p.cstate.cig.len();

    // then produce the sequence representing the insertion
    let mut k = p.cigar_index + 1;
    let mut offset = 1;
    while k < ncig {
        match p.cstate.cig[k] {
            Cigar::Pad(l) => seq_buf.extend(std::iter::repeat_n(b'*', l as usize)),
            Cigar::Ins(l) => {
                for _ in 0..l as usize {
                    read_pos = p.qpos + offset - p.del as usize;
                    read_base = p.rec.seq()[read_pos];
                    offset += 1;
                    seq_buf.push(read_base.to_ascii_uppercase());
                }
            }
            Cigar::Del(l) => {
                *ndel = l as i32;
                break;
            }
            _ => break,
        }
        k += 1;
    }

    Ok(())
}

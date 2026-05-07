use crate::errors::{Error, ErrorKind};
use crate::{alignment::PileupAlignment, output::OrderedPileupOutput, refseq::RefSeqHandle};
use indexmap::IndexMap;
use rust_htslib::bam::record::Cigar;
use std::ops::AddAssign;

#[derive(Clone, Debug)]
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
}

unsafe impl Sync for BaseDepthString {}
unsafe impl Send for BaseDepthString {}

impl OrderedPileupOutput for BaseDepthString {
    fn tid(&self) -> i32 {
        self.tid
    }

    fn pos(&self) -> i64 {
        self.pos
    }

    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, _ref_seq: &RefSeqHandle) {
        self.update(tid, pos, ref_name);
    }

    #[inline(always)]
    fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error> {
        self.intake(p, refseq)
    }

    #[inline(always)]
    fn write<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Error> {
        let mut buf = itoa::Buffer::new();

        writer.write_all(self.ref_name.as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.pos + 1).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.depth).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.a).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.g).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.c).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.t).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.n).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.gap).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.refskip).as_bytes())?;

        writer.write_all(b"\t[")?;

        let n = self.insertions.len() - 1;
        for (i, (ins, cnt)) in self.insertions.iter().enumerate() {
            writer.write_all(buf.format(*cnt).as_bytes())?;
            writer.write_all(ins)?;
            if i < n {
                writer.write_all(b" ")?
            };
        }

        writer.write_all(b"]\t[")?;

        let n = self.deletions.len() - 1;
        for (i, (del, cnt)) in self.deletions.iter().enumerate() {
            writer.write_all(buf.format(*cnt).as_bytes())?;
            writer.write_all(del)?;
            if i < n {
                writer.write_all(b" ")?
            };
        }

        write!(writer, "]")?;
        writeln!(writer)?;

        Ok(())
    }

    #[inline(always)]
    fn depth(&self) -> u32 {
        self.depth
    }

    #[inline(always)]
    fn clear(&mut self) {
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

    fn new() -> Self {
        Self::new()
    }
}

#[allow(clippy::new_without_default)]
impl BaseDepthString {
    pub fn new() -> Self {
        Self {
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

    pub fn update(&mut self, tid: i32, ref_pos: i64, ref_name: &str) {
        self.tid = tid;
        self.pos = ref_pos;
        if self.ref_name != ref_name {
            self.ref_name = ref_name.to_string()
        }
    }

    #[inline(always)]
    pub fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error> {
        self.depth += 1;
        self.register_pileup(p, refseq)
    }

    #[inline(always)]
    pub fn register_pileup(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error> {
        match p.del {
            false => {
                let readbase = if p.qpos < p.rec.seq_len() {
                    p.rec.seq()[p.qpos]
                } else {
                    b'n'
                };

                match readbase.to_ascii_uppercase() {
                    b'A' => self.a += 1,
                    b'G' => self.g += 1,
                    b'C' => self.c += 1,
                    b'T' => self.t += 1,
                    b'N' => self.n += 1,
                    other => {
                        return Err(Error::from(ErrorKind::AnomalousData(format!(
                            "Unrecognized nucleotide character: {other}"
                        ))))
                    }
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
                refbase = if let Some(refseq) = refseq.as_ref() {
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

#[inline(always)]
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
                    read_base = if read_pos < p.rec.seq_len() {
                        p.rec.seq()[read_pos]
                    } else {
                        b'n'
                    };
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

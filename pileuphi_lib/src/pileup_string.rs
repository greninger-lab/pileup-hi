use crate::alignment::PileupAlignment;
use crate::errors::Error;
use crate::output::OrderedPileupOutput;
use crate::refseq::RefSeqHandle;
use rust_htslib::bam::record::Cigar;
use std::io::Write;

const LAST_POS: u8 = b'$';
const FIRST_POS: u8 = b'^';

const F_MATCH: u8 = b'.';
const R_MATCH: u8 = b',';

const F_REFSKIP: u8 = b'>';
const R_REFSKIP: u8 = b'<';

#[derive(Clone, Debug)]
/// The standard samtools mpileup string. See the samtools mpileup docs for details.
pub struct PileupString {
    seqbuf: Vec<u8>,
    qualbuf: Vec<u8>,
    tid: i32,
    ref_name: String,
    ref_pos: i64,
    ref_base: u8,
    pub depth: u32,
}

unsafe impl Send for PileupString {}
unsafe impl Sync for PileupString {}

impl OrderedPileupOutput for PileupString {
    fn tid(&self) -> i32 {
        self.tid
    }

    fn pos(&self) -> i64 {
        self.ref_pos
    }

    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, ref_seq: &RefSeqHandle) {
        self.update(tid, pos, ref_name, ref_seq);
    }

    fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error> {
        self.intake(p, refseq)
    }

    fn write<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Error> {
        self.write(writer)
    }

    fn depth(&self) -> u32 {
        self.depth
    }

    fn new() -> Self {
        Self::new()
    }

    fn clear(&mut self) {
        self.depth = 0;
        self.seqbuf.clear();
        self.qualbuf.clear();
    }
}

impl PileupString {
    pub fn update(&mut self, tid: i32, ref_pos: i64, ref_name: &str, ref_seq: &RefSeqHandle) {
        self.tid = tid;
        self.ref_pos = ref_pos;

        self.ref_base = if let Some(seq) = ref_seq.as_ref() {
            *seq.get(ref_pos as usize).unwrap_or(&b'N')
        } else {
            b'N'
        };

        if self.ref_name != ref_name {
            self.ref_name = ref_name.to_string();
        }
    }

    #[inline(always)]
    pub fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error> {
        self.depth += 1;
        write_plp(p, self.ref_pos, &mut self.seqbuf, &mut self.qualbuf, refseq)?;
        Ok(())
    }

    #[inline(always)]
    pub fn write<W: std::io::Write>(&self, writer: &mut W) -> Result<(), Error> {
        let mut buf = itoa::Buffer::new();
        // print! {"{}\t{}\t{}\t{}\t", self.ref_name, self.ref_pos + 1, char::from(self.ref_base), self.depth }
        writer.write_all(self.ref_name.as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.ref_pos + 1).as_bytes())?;
        writer.write_all(b"\t")?;

        writer.write_all(&[self.ref_base])?;
        writer.write_all(b"\t")?;

        writer.write_all(buf.format(self.depth).as_bytes())?;
        writer.write_all(b"\t")?;

        if self.seqbuf.is_empty() {
            // write!(writer, "*\t")?
            writer.write_all(b"*\t")?
            // print! {"*\t"}
        } else {
            // write!(writer, "{}\t", std::str::from_utf8_unchecked(&self.seqbuf))?;
            writer.write_all(&self.seqbuf)?;
            writer.write_all(b"\t")?;
            // print! {"{}\t", std::str::from_utf8(&self.seqbuf)?}
        }

        if self.qualbuf.is_empty() {
            // write!(writer, "*")?
            writer.write_all(b"*")?;
        } else {
            // print! {"{}", std::str::from_utf8(&self.qualbuf)?}
            // unsafe { write!(writer, "{}", std::str::from_utf8_unchecked(&self.qualbuf))? }
            writer.write_all(&self.qualbuf)?;
        }

        writeln!(writer)?;

        Ok(())
    }

    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self {
            tid: 0,
            ref_pos: 0,
            ref_base: 0,
            depth: 0,
            ref_name: "".to_string(),
            qualbuf: Vec::new(),
            seqbuf: Vec::new(),
        }
    }
}

// cap qualitites at max of 126; this also helps avoid non-ascii output
pub fn get_qual(qual: u8) -> u8 {
    match qual.cmp(&93).is_gt() {
        true => 126,
        false => qual + 33,
    }
}

// TODO: take arguments that determine verbosity of reported insertion, e.g. full sequence or just
// length?
#[inline(always)]
pub fn expand_insertions(
    p: &PileupAlignment,
    seq_buf: &mut Vec<u8>,
    ndel: &mut i32,
    decorate: bool,
) -> Result<(), Error> {
    let is_rev = p.rec.is_reverse();
    let mut read_pos: usize;
    let mut read_base: u8;

    *ndel = 0;
    // first measure how many insertion bases there are.
    let mut len_indel = 0;
    let ncig = p.cstate.cig.len();
    let mut k = p.cigar_index + 1;
    while k < ncig {
        match p.cstate.cig[k] {
            Cigar::Pad(l) | Cigar::Ins(l) => {
                len_indel += l;
            }
            _ => break,
        }
        k += 1;
    }

    if decorate {
        write!(seq_buf, "+{}", len_indel)?
    };

    // then produce the sequence representing the insertion
    k = p.cigar_index + 1;
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
                    match is_rev {
                        true => seq_buf.push(read_base.to_ascii_lowercase()),
                        false => seq_buf.push(read_base.to_ascii_uppercase()),
                    }
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

#[inline(always)]
pub fn write_plp(
    p: &PileupAlignment,
    pos: i64,
    seq_buf: &mut Vec<u8>,
    qualbuf: &mut Vec<u8>,
    refseq: &RefSeqHandle,
) -> Result<(), Error> {
    if p.head {
        seq_buf.push(FIRST_POS);
        seq_buf.push(get_qual(p.rec.mapq()));
    }

    let is_rev = p.rec.is_reverse();
    let mut refbase: u8;

    match p.del {
        false => {
            refbase = if let Some(refseq) = refseq.as_ref() {
                refseq[pos as usize]
            } else {
                b'('
            };
            let readbase = if p.qpos < p.rec.seq_len() {
                p.rec.seq()[p.qpos]
            } else {
                b'n'
            };

            if readbase.eq_ignore_ascii_case(&refbase) {
                match is_rev {
                    true => seq_buf.push(R_MATCH),
                    false => seq_buf.push(F_MATCH),
                }
            } else {
                match is_rev {
                    true => seq_buf.push(readbase.to_ascii_lowercase()),
                    false => seq_buf.push(readbase.to_ascii_uppercase()),
                }
            }
        }

        true => {
            if p.refskip {
                match is_rev {
                    true => seq_buf.push(R_REFSKIP),
                    false => seq_buf.push(F_REFSKIP),
                };
            } else {
                seq_buf.push(b'*');
            }
        }
    }

    let mut del_len = -p.indel;
    if p.indel > 0 {
        expand_insertions(p, seq_buf, &mut del_len, true)?;
    }

    if del_len > 0 {
        write!(seq_buf, "{}", -del_len)?;
        for i in 1..=del_len as i64 {
            refbase = if let Some(refseq) = refseq.as_ref() {
                refseq[(pos + i) as usize]
            } else {
                b'N'
            };

            match is_rev {
                false => seq_buf.push(refbase.to_ascii_uppercase()),
                true => seq_buf.push(refbase.to_ascii_lowercase()),
            }
        }
    }

    if p.tail {
        seq_buf.push(LAST_POS);
    }

    // finally,we add PHRED qual
    qualbuf.push(get_qual(*p.rec.qual().get(p.qpos).unwrap_or(&0)));

    Ok(())
}

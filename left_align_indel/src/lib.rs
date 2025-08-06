use rust_htslib::bam::{
    Format, Header, Read, Reader, Record, Writer,
    ext::BamRecordExtensions,
    record::{Cigar, CigarString},
};

use bio::io::fasta::IndexedReader;

use std::path::Path;

/// Shift insertions and deletions spanning repetitive sequence to the left-most
/// position at which the reference sequence and read sequence agree (on the repetitive sequence).
///
/// This is done to mitigate arbitrary placement of cigars along a stretch of repeat sequence performed by some aligners.
/// Ideally, more consistent positioning of indels can lead to more accurate variant calling.
///
/// This is largely taken from GATK's LeftAlignIndels.java.
///
///
/// EXAMPLE:
/// ```text
/// REFERENCE: AGCTATATAGC
/// QUERY:     AGCTATA--GC
/// CIGAR: 7M2D2M
/// ```
/// we want to shift the deletion here to the very left of the read's and ref's repeat sequence.
/// So we can get this:
/// ```text
/// REFERENCE: AGCTATATAGC
/// QUERY:     AGC--TATAGC
/// CIGAR: 3M2D6M
/// ```
pub fn left_align_indels(record: &Record, refseq: &[u8]) -> CigarString {
    let seqlen = record.seq_len_from_cigar(false);
    let reflen = record.reference_end() - record.pos();

    let mut read_range = CoordinateRange {
        start: seqlen as i32,
        end: seqlen as i32,
    };

    let mut ref_range = CoordinateRange {
        start: (record.pos() as i32 + reflen as i32),
        end: (record.pos() as i32 + reflen as i32),
    };

    let mut cigar_reverse: Vec<Cigar> = Vec::with_capacity(record.cigar_len());
    let mut op: Cigar;
    let read_seq = record.seq().as_bytes();

    let n = record.cigar_len();
    for i in (0..n).rev() {
        op = record.cigar()[i as usize];

        if is_indel(op) {
            // deletion found, so expand boundary to hold indel bases
            read_range.shift_start_left(length_on_read(op));
            ref_range.shift_start_left(length_on_ref(op));
        } else if read_range.len() + ref_range.len() == 0 {
            // no indel found yet, so just shift op by left and record current cigar op
            cigar_reverse.push(op);
            read_range.shift_left(length_on_read(op));
            ref_range.shift_left(length_on_ref(op));
        } else {
            // we've hit a non-indel to the left of an indel. Attempt to shift left.
            let maxshift = if is_alignment(op) { op.len() as i32 } else { 0 };
            let (sequences, bounds) = (
                vec![refseq, read_seq.as_slice()],
                vec![&mut ref_range, &mut read_range],
            );
            let (lshift, rshift) = shift_indel_left(sequences, bounds, maxshift, true);

            // add new match sequences, if we actually shifted the deletion to the left in
            // fuse with previous match if it exists
            if let Some(lastop) = cigar_reverse.last_mut() {
                match lastop {
                    Cigar::Match(l) => *lastop = Cigar::Match(*l + rshift as u32),
                    _ => (),
                };
            } else {
                cigar_reverse.push(Cigar::Match(rshift as u32));
            }

            let emit_indel = i == 0 || lshift < maxshift || !is_alignment(op);

            // the lshift to use for later
            // make sure i32 conversion is okay here
            // we shifted to the right in some cases.
            let left_new_match = if lshift < 0 { -1 * lshift } else { 0 };

            // if we didn't shift left at all, all of the remaining match op is still left.
            let remaining_left_bases = if lshift < 0 {
                op.len() as i32
            } else {
                op.len() as i32 - lshift
            };

            if emit_indel {
                if ref_range.len() > 0 {
                    cigar_reverse.push(Cigar::Del(ref_range.len() as u32));
                }
                if read_range.len() > 0 {
                    cigar_reverse.push(Cigar::Ins(read_range.len() as u32));
                }

                // advance the interval by the len of the op we just emitted
                ref_range.shift_end_left(ref_range.len());
                read_range.shift_end_left(read_range.len());

                // advance the coordinate range by the remainder of the current match op
                ref_range.shift_left(
                    left_new_match
                        + length_on_ref(op)
                            .eq(&0)
                            .then_some(0)
                            .unwrap_or(remaining_left_bases),
                );
                read_range.shift_left(
                    left_new_match
                        + length_on_read(op)
                            .eq(&0)
                            .then_some(0)
                            .unwrap_or(remaining_left_bases),
                );
            }

            // println! {"shifts: {lshift} {rshift}"}
            // println! {"bases remaining on left: {remaining_left_bases}"}
            // println! {"new match on left due to trimming: {left_new_match}"}

            // now add the match op to the cigar string
            // note that we aren't preserving the exact cigar (e.g. X, Eq, Match), just using
            // match for now.
            if left_new_match + remaining_left_bases > 0 {
                cigar_reverse.push(Cigar::Match(
                    left_new_match as u32 + remaining_left_bases as u32,
                ));
            }
        }
    }

    assert_eq!(
        read_range.start,
        0,
        "{:?}",
        std::str::from_utf8(record.qname()).unwrap()
    );
    cigar_reverse.reverse();
    return CigarString::from(cigar_reverse);
}

#[derive(Debug)]
pub struct CoordinateRange {
    start: i32,
    end: i32,
}

impl CoordinateRange {
    pub fn len(&self) -> i32 {
        self.end - self.start
    }

    pub fn shift_start_left(&mut self, shift: i32) {
        self.start -= shift;
    }

    pub fn shift_start_right(&mut self, shift: i32) {
        self.start += shift;
    }

    pub fn shift_end_left(&mut self, shift: i32) {
        self.end -= shift;
    }

    pub fn shift_end_right(&mut self, shift: i32) {
        self.end += shift;
    }

    pub fn shift_left(&mut self, shift: i32) {
        self.start -= shift;
        self.end -= shift;
    }

    pub fn shift_right(&mut self, shift: i32) {
        self.start += shift;
        self.end += shift;
    }
}

/* Given a set of intervals denoting an indel on the given sequences:
* 1. converge on indel end in all intervals, trimming off matching bases on either end if trim
*    is enabled.
* 2. perform a sliding window along the reference and read sequence, shifting the coordinates of
*    the selected indel to the left as long as its sequence is repeated 5' in all sequences
*    supplied.
*
* 3. return the amount by which the left and right bounds of all interval have been shifted.
*/
pub fn shift_indel_left(
    sequences: Vec<&[u8]>,
    mut bounds: Vec<&mut CoordinateRange>,
    maxshift: i32,
    trim: bool,
) -> (i32, i32) {
    let mut startshift: i32 = 0;
    let mut endshift: i32 = 0;

    let mut minsize = bounds.iter().map(|x| x.len()).min().unwrap() as usize;

    while trim && minsize > 0 && last_base_on_right_same(&sequences, &bounds) {
        bounds.iter_mut().for_each(|b| b.shift_end_left(1));
        minsize -= 1;
        endshift += 1;
    }

    while trim && minsize > 0 && first_base_on_left_same(&sequences, &bounds) {
        bounds.iter_mut().for_each(|b| b.shift_start_left(1));
        minsize -= 1;
        startshift -= 1;
    }

    while startshift < maxshift
        && next_base_on_left_same(&sequences, &bounds)
        && last_base_on_right_same(&sequences, &bounds)
    {
        bounds.iter_mut().for_each(|b| b.shift_left(1));
        startshift += 1;
        endshift += 1;
    }

    (startshift, endshift)
}

pub fn last_base_on_right_same(sequences: &Vec<&[u8]>, bounds: &Vec<&mut CoordinateRange>) -> bool {
    let last_right_base = sequences[0][bounds[0].end as usize - 1];

    for i in 0..sequences.len() {
        if sequences[i][bounds[i].end as usize - 1] != last_right_base {
            return false;
        }
    }
    true
}

pub fn first_base_on_left_same(sequences: &Vec<&[u8]>, bounds: &Vec<&mut CoordinateRange>) -> bool {
    let first_left_base = sequences[0][bounds[0].start as usize];
    for i in 0..sequences.len() {
        if sequences[i][bounds[i].start as usize] != first_left_base {
            return false;
        }
    }

    true
}

pub fn next_base_on_left_same(sequences: &Vec<&[u8]>, bounds: &Vec<&mut CoordinateRange>) -> bool {
    let next_base_left = sequences[0][bounds[0].start as usize - 1];

    for i in 0..sequences.len() {
        if sequences[i][bounds[i].start as usize - 1] != next_base_left {
            return false;
        }
    }
    true
}

fn is_indel(cig: Cigar) -> bool {
    match cig {
        Cigar::Ins(_) | Cigar::Del(_) => true,
        _ => false,
    }
}

fn is_alignment(cig: Cigar) -> bool {
    match cig {
        Cigar::Equal(_) | Cigar::Match(_) | Cigar::Diff(_) => true,
        _ => false,
    }
}

fn length_on_ref(cig: Cigar) -> i32 {
    if is_alignment(cig) || matches!(cig, Cigar::Del(_)) {
        cig.len() as i32
    } else {
        0
    }
}

fn length_on_read(cig: Cigar) -> i32 {
    if is_alignment(cig) || matches!(cig, Cigar::Ins(_)) || matches!(cig, Cigar::SoftClip(_)) {
        cig.len() as i32
    } else {
        0
    }
}

pub fn bam_from_file(file: &str) -> Reader {
    let mut reader = rust_htslib::bam::Reader::from_path(file).unwrap();
    reader.set_threads(8).unwrap();

    reader
}

pub fn make_bam_out(file: &str, reader: &Reader) -> Writer {
    let mut writer =
        Writer::from_path(file, &Header::from_template(&reader.header()), Format::Sam).unwrap();
    writer.set_threads(8).unwrap();

    writer
}

pub fn load_seq(fasta: &str, seqid: Option<&str>) -> Vec<u8> {
    let mut out = vec![];
    let mut reader = IndexedReader::from_file(&Path::new(fasta)).unwrap();
    if let Some(seqid) = seqid {
        reader.fetch_all(seqid).unwrap();
    } else {
        reader.fetch_all_by_rid(0).unwrap();
    }
    reader.read(&mut out).unwrap();

    return out;
}

#[cfg(test)]
mod tests {
    use super::*;
    use rust_htslib::bam::{
        Read, Record,
        ext::BamRecordExtensions,
        record::{Cigar, CigarString, CigarStringView},
    };

    use std::{
        fs::File,
        io::{BufRead, BufReader},
    };

    pub fn assert_files_equal(f1: &str, f2: &str, comment: &str) {
        let lines1 = BufReader::new(File::open(f1).unwrap()).lines();
        let lines2 = BufReader::new(File::open(f2).unwrap()).lines();
        let mut l1 @ mut l2: String;

        for (_l1, _l2) in lines1.zip(lines2) {
            l1 = _l1.unwrap();
            l2 = _l2.unwrap();

            if l1.starts_with(comment) {
                continue;
            }

            assert_eq!(l1, l2, "{} {}", l1, l2);
        }
    }

    #[test]
    fn test_contrived() {
        let reff = b"AGCTATATAGC";
        let quer = b"AGCTATAGC";
        let qual = b"#########";

        let cig = CigarString::from(vec![Cigar::Match(7), Cigar::Del(2), Cigar::Match(2)]);

        let mut r = Record::new();
        r.set(b"1", Some(&cig), quer, qual);
        r.set_pos(0);
        r.unset_unmapped();

        let seqlen = r.seq_len_from_cigar(false);
        let reflen = r.reference_end() - r.pos();

        assert_eq!(seqlen, 9);
        assert_eq!(reflen, 11);

        let cig = left_align_indels(&r, reff);

        assert_eq!(
            cig,
            CigarString::from(vec![Cigar::Match(3), Cigar::Del(2), Cigar::Match(6)])
        );

        println! {"{cig}"}
    }

    #[test]
    fn test_real() {
        let mut input = bam_from_file("non-aligned.sam");
        let mut output = make_bam_out("realigned.sam", &input);

        println! {"getting refseq"}
        let refseq = load_seq("Homo_sapiens_assembly19.fasta", Some("20"));
        let mut old_cig: CigarStringView;
        let mut cig: CigarString;

        for mut r in input.records().flatten() {
            old_cig = r.cigar();
            if r.is_unmapped() {
                continue;
            }
            assert_eq!(r.tid(), 0);
            // println!{"{:?} {}", r.seq(), r.cigar()}
            cig = left_align_indels(&r, &refseq);
            println! {"{old_cig} | {cig}"}
            r.set_cigar(Some(&cig));
            output.write(&r).unwrap();
        }
    }
}

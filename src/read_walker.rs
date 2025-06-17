use rust_htslib::bam::record::Cigar;
use rust_htslib::bam::Record;

pub struct IterCigarMatches {
    cigar: Vec<Cigar>,
    cigar_index: usize,
    match_remaining: u32,
    read_pos: usize,
    pub genome_pos: i64,
    found: bool,
}

pub trait WalkMatches {
    fn walk_matches(&self) -> IterCigarMatches;
}

impl WalkMatches for Record {
    fn walk_matches(&self) -> IterCigarMatches {
        IterCigarMatches {
            cigar: self.cigar().take().0,
            cigar_index: 0,
            match_remaining: 0,
            read_pos: 0,
            genome_pos: self.pos(),
            found: false,
        }
    }
}

impl Iterator for IterCigarMatches {
    type Item = (usize, i64);
    fn next(&mut self) -> Option<Self::Item> {
        if self.match_remaining > 0 {
            self.match_remaining -= 1;
            self.genome_pos += 1;
            self.read_pos += 1;
            return Some((self.read_pos - 1, self.genome_pos - 1));
        }

        loop {
            if self.cigar_index >= self.cigar.len() {
                // self.cigar_index -= 1; // don't use iterator after this
                return None;
            }

            let entry = self.cigar[self.cigar_index];

            match entry {
                Cigar::Equal(l) | Cigar::Match(l) | Cigar::Diff(l) => {
                    self.found = true;
                    self.read_pos += 1;
                    self.genome_pos += 1;
                    self.match_remaining = l - 1;
                    self.cigar_index += 1;
                    return Some((self.read_pos - 1, self.genome_pos - 1));
                }
                Cigar::Ins(l) | Cigar::SoftClip(l) => {
                    self.read_pos += l as usize;
                    self.cigar_index += 1;
                }
                Cigar::Del(l) => {
                    self.genome_pos += l as i64;
                    self.cigar_index += 1;
                }

                Cigar::RefSkip(_) | Cigar::Pad(_) | Cigar::HardClip(_) => (),
            }
        }
    }
}

impl IterCigarMatches {
    pub fn after_del(&self) -> bool {
        if self.cigar_index <= 1 {
            false
        } else {
            match self.cigar[self.cigar_index - 2] {
                Cigar::Del(_) => true,
                _ => false,
            }
        }
    }

    pub fn _cur_cigar(&self) -> Cigar {
        self.cigar[self.cigar_index]
    }
}

#[test]
fn test_walk1() {
    let mut r = Record::new();
    r.set(
        b"r1",
        Some(&CigarString(vec![Cigar::Match(5)])),
        b"GATGA",
        b"#####",
    );

    r.set_pos(0);

    let ret = r.walk_matches().collect::<Vec<(usize, i64)>>();
    assert_eq!(ret, vec![(0, 0), (1, 1), (2, 2), (3, 3), (4, 4)])
}

#[test]
fn test_walk2() {
    let mut r = Record::new();
    r.set(
        b"r2",
        Some(&CigarString(vec![
            Cigar::Match(4),
            Cigar::Del(1),
            Cigar::Match(1),
        ])),
        b"GATGA",
        b"#####",
    );

    r.set_pos(0);

    let mut walker = r.walk_matches();

    let mut ret = vec![];

    while let Some(next) = walker.next() {
        ret.push(next);
    }

    assert_eq!(ret, vec![(0, 0), (1, 1), (2, 2), (3, 3), (4, 5)]);

    assert!(matches!(
        walker.cigar[walker.cigar_index - 2],
        Cigar::Del(1)
    ));

    println! {"{}", walker.cigar[walker.cigar_index - 2]}
    assert!(walker.after_del());
}

#[cfg(test)]
use rust_htslib::bam::record::CigarString;

#[test]
fn test_walk3() {
    let mut r = Record::new();
    r.set(
        b"r2",
        Some(&CigarString(vec![
            Cigar::Ins(1),
            Cigar::Match(4),
            Cigar::Del(3),
            Cigar::Match(2),
        ])),
        b"GGATGAA",
        b"#######",
    );

    // -012345678 | ref pos
    // *GATG---AA
    // 0123456789 | read pos

    r.set_pos(0);

    let mut walker = r.walk_matches();

    let mut ret = vec![];

    while let Some(next) = walker.next() {
        ret.push(next);
    }

    assert_eq!(ret, vec![(1, 0), (2, 1), (3, 2), (4, 3), (5, 7), (6, 8)]);

    assert!(matches!(
        walker.cigar[walker.cigar_index - 2],
        Cigar::Del(3)
    ));

    assert!(walker.after_del());

    println! {"{}", walker.cigar[walker.cigar_index - 2]}
}

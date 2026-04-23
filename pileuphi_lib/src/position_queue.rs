use crate::errors::{Error, ErrorKind};
use log::warn;
use rust_htslib::bam::HeaderView;
use std::ops::Not;

/// A raw pileup region not yet validated to actually exist in a BAM header.
pub struct RawPileupRegion {
    name: String,
    start: i64,
    end: i64,
}

/// A pileup region with associated header TID, presumably validated.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct GenomeInterval {
    pub tid: i64,
    pub name: String,
    pub start: i64,
    pub end: i64,
}

pub struct GenomeIntervalIterator<'a> {
    cur_start: i64,
    cur_end: i64,
    chunk_size: i64,
    interval: &'a GenomeInterval,
    name: String,
    exhausted: bool,
}

impl<'a> GenomeIntervalIterator<'a> {
    pub fn new(interval: &'a GenomeInterval, chunk_size: i64) -> Self {
        Self {
            cur_start: interval.start,
            cur_end: interval.start + chunk_size,
            chunk_size,
            interval,
            name: interval.name.to_string(),
            exhausted: false,
        }
    }
}

impl Iterator for GenomeIntervalIterator<'_> {
    type Item = GenomeInterval;

    fn next(&mut self) -> Option<Self::Item> {
        let ret: Option<Self::Item>;
        if self.exhausted {
            return None;
        };

        if self.chunk_size > self.interval.end - self.interval.start {
            self.exhausted = true;
            return Some(self.interval.clone());
        }

        // start and end of next chunk is still within original interval
        if self.cur_end < self.interval.end {
            ret = Some(GenomeInterval {
                tid: self.interval.tid,
                name: self.name.to_string(),
                start: self.cur_start,
                end: self.cur_end,
            });

        // start is still within but end is outside. Clamp end to max coordinate.
        } else if self.cur_start < self.interval.end {
            ret = Some(GenomeInterval {
                tid: self.interval.tid,
                name: self.name.to_string(),
                start: self.cur_start,
                end: std::cmp::min(self.cur_end, self.interval.end),
            });

        // completely outside window
        } else {
            self.exhausted = true;
            ret = None
        }

        // advance and yield
        self.cur_start += self.chunk_size;
        self.cur_end += self.chunk_size;

        ret
    }
}

impl GenomeInterval {
    #[allow(dead_code)]
    pub fn chunks(&self, chunk_size: i64) -> GenomeIntervalIterator<'_> {
        GenomeIntervalIterator::new(self, chunk_size)
    }

    #[allow(dead_code)]
    pub fn n_chunks(&self, n_chunks: i64) -> GenomeIntervalIterator<'_> {
        if n_chunks == 1 {
            GenomeIntervalIterator::new(self, i64::MAX)
        } else {
            GenomeIntervalIterator::new(self, (self.end - self.start) / n_chunks + 1)
        }
    }

    pub fn is_empty(&self) -> bool {
        self.end >= self.start
    }

    pub fn len(&self) -> i64 {
        self.end - self.start
    }
}

// split a str by a delimiter and convert empty prefix/suffix to None
fn split_check_ends(s: &str, delim: char) -> Option<(Option<&str>, Option<&str>)> {
    s.split_once(delim).map(|(pre, post)| {
        (
            pre.is_empty().not().then_some(pre),
            post.is_empty().not().then_some(post),
        )
    })
}

/// Parse any string for being compliant for the SAM region format, e.g.
/// chr1:400-801
fn parse_region_string(s: &str) -> Result<RawPileupRegion, Error> {
    match split_check_ends(s, ':') {
        // no coordinates, just entire reference
        None => Ok(RawPileupRegion {
            name: s.to_string(),
            start: 0,
            end: i64::MAX,
        }),
        Some((None, _)) => Err(Error::from(ErrorKind::BadInputRegions(format!(
            "Invalid region string {s}: ref name must come before ':'"
        )))),

        Some((Some(_), None)) => Err(Error::from(ErrorKind::BadInputRegions(format!(
            "Invalid region string {s}: coordinates must come after ':'"
        )))),

        // parse coordinates
        Some((Some(ref_str), Some(pos_str))) => match split_check_ends(pos_str, '-') {
            None => {
                // no dashes, meaning we expect a one-coordinate interval: e.g. Chr1:400
                let pos = pos_str.replace(",", "").parse::<i64>()?;

                if pos < 0 {
                    return Err(Error::from(ErrorKind::BadInputRegions(format!(
                        "Cannot have negative pos {pos}: {s}"
                    ))));
                };

                Ok(RawPileupRegion {
                    name: ref_str.to_string(),
                    start: pos - 1,
                    end: i64::MAX,
                })
            }

            Some((None, _)) => Err(Error::from(ErrorKind::BadInputRegions(format!(
                "Invalid region string {s}: must have start coordinate before '-'"
            )))),

            Some((Some(_), None)) => Err(Error::from(ErrorKind::BadInputRegions(format!(
                "Invalid region string {s}: must have end coordinate after '-'"
            )))),

            // We have something before and after a dash, so we expect numbers on both sides...
            Some((Some(start), Some(end))) => {
                let mut start_pos = start.replace(",", "").parse::<i64>()?;
                let end_pos = end.replace(",", "").parse::<i64>()?;

                if start_pos < 0 {
                    return Err(Error::from(ErrorKind::BadInputRegions(format!(
                        "Invalid region string {s}: cannot have negative start pos: {start_pos}"
                    ))));
                };

                if start_pos == 0 {
                    warn!("start pos is 0: input regions are 1-based. Setting start pos to 1");
                    start_pos = 1;
                }

                if end_pos < 0 {
                    return Err(Error::from(ErrorKind::BadInputRegions(format!(
                        "Invalid region string {s}: cannot have negative end pos: {end_pos}"
                    ))));
                };
                if end_pos < start_pos {
                    return Err(Error::from(ErrorKind::BadInputRegions(format!(
                        "Invalid region string {s}: cannot have start end pos ({end_pos}) be smaller than start pos ({start_pos})"
                    ))));
                };

                Ok(RawPileupRegion {
                    name: ref_str.to_string(),
                    start: start_pos - 1,
                    end: end_pos - 1,
                })
            }
        },
    }
}

pub fn create_region_queue(argstr: &str, header: &HeaderView) -> Result<Vec<GenomeInterval>, Error> {
    let mut rawregions = vec![];

    for s in argstr.split_terminator(",") {
        rawregions.push(parse_region_string(s)?);
    }

    intervals_from_regions(header, rawregions)
}

/// Create a PositionQueue from the entire header
pub fn intervals_from_header(header: &HeaderView) -> Result<Vec<GenomeInterval>, Error> {
    let mut queue = Vec::new();

    for tid in 0..header.target_count() {
        let end = header
            .target_len(tid)
            .ok_or(Error::from(ErrorKind::AnomalousData(format!(
                "Faulty BAM header: unable to find target len for ref ID {tid}"
            ))))? as i64;

        let name = header.tid2name(tid);

        let reg = GenomeInterval {
            name: std::str::from_utf8(name)?.to_string(),

            tid: i64::from(tid),
            start: 0,
            end,
        };

        queue.push(reg);
    }

    Ok(queue)
}

/// Create a position queue from a list of pileup regions, validating to make sure they agree
/// with the given SAM header.
pub fn intervals_from_regions(
    header: &HeaderView,
    regions: Vec<RawPileupRegion>,
) -> Result<Vec<GenomeInterval>, Error> {
    if regions.is_empty() {
        return Err(Error::from(ErrorKind::AnomalousData(
            "List of BAM regions is empty!".to_string(),
        )));
    }

    let tnames: Vec<&str> = header
        .target_names()
        .into_iter()
        .map(|s| std::str::from_utf8(s))
        .collect::<Result<Vec<&str>, _>>()?;

    let mut queue = Vec::new();

    for rawreg in regions {
        let mut found = false;
        for (tid, canonname) in tnames.iter().enumerate() {
            if rawreg.name == *canonname {
                let canonlen = header
                    .target_len(tid as u32)
                    .ok_or(Error::from(ErrorKind::AnomalousData(format!(
                        "Faulty BAM header: unable to find target len for ref ID {tid}"
                    ))))?;

                let end = if rawreg.end > canonlen as i64 {
                    canonlen as i64
                } else {
                    rawreg.end + 1
                };

                found = true;
                queue.push(GenomeInterval {
                    name: rawreg.name.clone(),
                    tid: i64::try_from(tid)?,
                    start: rawreg.start,
                    end,
                    // end: rawreg.end.min(canonlen as i64),
                })
            }
        }

        if !found {
            return Err(Error::from(ErrorKind::AnomalousData(format!(
                "Couldn't find region name {} in BAM header",
                rawreg.name
            ))));
        }
    }

    Ok(queue)
}

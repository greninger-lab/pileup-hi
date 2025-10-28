use anyhow::{Context, Error};
use rust_htslib::bam::HeaderView;

/// A raw pileup region not yet validated to actually exist in a BAM header.
pub struct RawPileupRegion {
    name: String,
    start: i64,
    end: i64,
}

/// A pileup region with associated header TID, presumably validated.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct GenomeInterval {
    // pub name: String,
    pub tid: i64,
    pub start: i64,
    pub end: i64,
}

pub struct GenomeIntervalIterator<'a> {
    cur_start: i64,
    cur_end: i64,
    chunk_size: i64,
    interval: &'a GenomeInterval,
}

impl<'a> GenomeIntervalIterator<'a> {
    pub fn new(interval: &'a GenomeInterval, chunk_size: i64) -> Self {
        Self {
            cur_start: interval.start,
            cur_end: interval.start + chunk_size,
            chunk_size,
            interval,
        }
    }
}

impl Iterator for GenomeIntervalIterator<'_> {
    type Item = GenomeInterval;

    fn next(&mut self) -> Option<Self::Item> {
        let ret: Option<Self::Item>;

        // start and end of next chunk is still within original interval
        if self.cur_end < self.interval.end {
            ret = Some(GenomeInterval {
                tid: self.interval.tid,
                start: self.cur_start,
                end: self.cur_end,
            });

        // start is still within but end is outside. Clamp end to max coordinate.
        } else if self.cur_start < self.interval.end {
            ret = Some(GenomeInterval {
                tid: self.interval.tid,
                start: self.cur_start,
                end: std::cmp::min(self.cur_end, self.interval.end),
            });

        // completely outside window
        } else {
            ret = None
        }

        // advance and yield
        self.cur_start += self.chunk_size;
        self.cur_end += self.chunk_size;
        ret
    }
}

impl GenomeInterval {
    pub fn chunks(&self, chunk_size: i64) -> GenomeIntervalIterator<'_> {
        GenomeIntervalIterator::new(self, chunk_size)
    }
}

/// Parse any string for being compliant for the SAM region format, e.g.
/// chr1:400-801
fn parse_region_string(s: &str) -> Result<RawPileupRegion, Error> {
    // region strings should have one colon and a single dash
    let col_count = s.chars().filter(|c| *c == ':').count();
    if col_count != 1 {
        anyhow::bail!("Invalid number of colons ({col_count}) in region string");
    }

    let dash_count = s.chars().filter(|c| *c == '-').count();
    if dash_count != 1 {
        anyhow::bail!("Invalid number of dashes ({dash_count}) in region string");
    }

    let (ref_name, pos_str) = s.split_once(":").unwrap();
    let (start, end) = pos_str.split_once("-").unwrap();

    let mut start = start
        .parse::<i64>()
        .with_context(|| format!("Non-numeric start position: {}", start))?;

    // subtract because BAM positions start at zero.
    start -= 1;

    if start < 0 {
        anyhow::bail!("Cannot have a negative start to a region: {}", start);
    }

    let end = end
        .parse::<i64>()
        .with_context(|| format!("Non-numeric end position: {}", end))?;

    Ok(RawPileupRegion {
        name: ref_name.to_string(),
        start,
        end,
    })
}

fn parse_region_arg(s: &str) -> Result<Vec<RawPileupRegion>, Error> {
    s.split_terminator(' ').map(parse_region_string).collect()
}

pub fn create_region_queue(argstr: &str, header: &HeaderView) -> Result<PositionQueue, Error> {
    let rawregions = parse_region_arg(argstr)?;
    PositionQueue::new_from_regions(header, rawregions)
}

pub struct PositionQueue {
    pub queue: Vec<GenomeInterval>,
}

impl PositionQueue {
    /// Create a PositionQueue from the entire header
    pub fn new(header: &HeaderView) -> Result<Self, Error> {
        let mut queue = Vec::new();

        for tid in 0..header.target_count() {
            let end = header
                .target_len(tid)
                .context("Unable to get target len")?
                .try_into()?;

            let reg = GenomeInterval {
                // name: name.to_string(),
                tid: i64::from(tid),
                start: 0,
                end,
            };

            queue.push(reg);
        }

        Ok(Self { queue })
    }

    /// Create a position queue from a list of pileup regions, validating to make sure they agree
    /// with the given SAM header.
    pub fn new_from_regions(
        header: &HeaderView,
        regions: Vec<RawPileupRegion>,
    ) -> Result<Self, Error> {
        if regions.is_empty() {
            anyhow::bail!("Cannot supply empty regions list to PositionQueue builder!");
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
                        .target_len(u32::try_from(tid)?)
                        .context("Unable to get ref len")?;

                    if rawreg.end >= i64::try_from(canonlen)? {
                        anyhow::bail!(
                            "Supplied region end exceeds length of reference in header: {} vs {}",
                            rawreg.end,
                            canonlen
                        );
                    }

                    found = true;
                    queue.push(GenomeInterval {
                        // name: rawreg.name.clone(),
                        tid: i64::try_from(tid)?,
                        start: rawreg.start,
                        end: rawreg.end,
                    })
                }
            }

            if !found {
                anyhow::bail!("Unable to find reference {} in header!", rawreg.name);
            }
        }

        Ok(Self { queue })
    }
}

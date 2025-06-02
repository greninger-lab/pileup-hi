use anyhow::Error;
use rust_htslib::bam::{ext::BamRecordExtensions, Record};

const BAM_FPAIRED: u16 = 1;
const BAM_FPROPER_PAIR: u16 = 2;
const BAM_FMUNMAP: u16 = 8;
const BAM_FREVERSE: u16 = 16;
const BAM_FREAD1: u16 = 64;
const BAM_FREAD2: u16 = 128;
const BAM_FSECONDARY: u16 = 256;
const BAM_FQCFAIL: u16 = 512;
const BAM_FDUP: u16 = 1024;
const BAM_FSUPPLEMENTARY: u16 = 2048;

pub struct ReadFilter {
    inc_flag: u16,
    exc_flag: u16,
    count_orphans: bool,
    min_mapq: u8,
}

fn add_to_flag(flags: Vec<&str>, flag: &mut u16) -> Result<(), Error> {
    for f in flags.into_iter() {
        match f {
            "BAM_FPAIRED" => *flag |= BAM_FPAIRED,
            "BAM_FSECONDARY" => *flag |= BAM_FSECONDARY,
            "BAM_FPROPER_PAIR" => *flag |= BAM_FPROPER_PAIR,
            "BAM_FMUNMAP" => *flag |= BAM_FMUNMAP,
            "BAM_FREVERSE" => *flag |= BAM_FREVERSE,
            "BAM_FREAD1" => *flag |= BAM_FREAD1,
            "BAM_FREAD2" => *flag |= BAM_FREAD2,
            "BAM_SECONDARY" => *flag |= BAM_FSECONDARY,
            "BAM_FQCFAIL" => *flag |= BAM_FQCFAIL,
            "BAM_FDUP" => *flag |= BAM_FDUP,
            "BAM_FSUPPLEMENTARY" => *flag |= BAM_FSUPPLEMENTARY,
            _ => anyhow::bail!("Unrecognized BAM flag specified: {f}"),
        }
    }

    Ok(())
}

impl ReadFilter {
    pub fn new(min_mapq: u8, excl_flags: Vec<&str>, incl_flags: Vec<&str>) -> Result<Self, Error> {
        let mut s = Self {
            inc_flag: 0,
            exc_flag: 0,
            count_orphans: false,
            min_mapq: 0,
        };

        s.add_incl_flags(incl_flags)?;
        s.add_excl_flags(excl_flags)?;
        s.min_mapq = min_mapq;
        Ok(s)
    }
    pub fn add_excl_flags(&mut self, flags: Vec<&str>) -> Result<(), Error> {
        let mut flag: u16 = 0;
        add_to_flag(flags, &mut flag)?;
        self.exc_flag = flag;

        Ok(())
    }

    pub fn add_incl_flags(&mut self, flags: Vec<&str>) -> Result<(), Error> {
        let mut flag: u16 = 0;
        add_to_flag(flags, &mut flag)?;
        self.inc_flag = flag;

        Ok(())
    }

    pub fn check_read(&mut self, read: &Record) -> bool {
        let mut pass;

        // check if orphan pair
        if read.is_paired() && !read.is_proper_pair() && !self.count_orphans {
            return false;
        }

        if read.mapq() < self.min_mapq + 33 {
            return false;
        }

        // println! {"{} {:?} {:?}", read.inner.core.flag & self.inc_flag, read.inner.core.flag.to_ne_bytes(), self.inc_flag.to_ne_bytes()}

        // does the flag contain any bits that we filter for?
        pass = (read.inner.core.flag & self.inc_flag) < read.inner.core.flag;

        if !pass {
            return pass;
        }

        // println! {"{} {:?} {:?}", read.inner.core.flag & self.exc_flag, read.inner.core.flag.to_ne_bytes(), self.exc_flag.to_ne_bytes()}

        // does the flag contain any bits that we filter against?
        pass |= (read.inner.core.flag & self.exc_flag) == 0;

        pass
    }
}

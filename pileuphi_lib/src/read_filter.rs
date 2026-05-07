use crate::errors::Error;
use rust_htslib::bam::Record;

pub struct ReadFilter {
    inc_flag: u16,
    exc_flag: u16,
    count_orphans: bool,
}

#[derive(Clone, Copy)]
#[cfg(feature = "cli")]
#[derive(clap::ValueEnum)]
pub enum BamFlag {
    Paired = 1,
    ProperPair = 2,
    Unmapped = 4,
    MateUnmapped = 8,
    Reverse = 16,
    Read1 = 64,
    Read2 = 128,
    Secondary = 256,
    QCFail = 512,
    Duplicate = 1024,
    Supplementary = 2048,
}

// don't change these, since they need to match with the clap value enum strings
impl BamFlag {
    pub fn as_str(&self) -> &str {
        match self {
            BamFlag::Paired => "paired",
            BamFlag::ProperPair => "proper-pair",
            BamFlag::Unmapped => "unmapped",
            BamFlag::MateUnmapped => "mate-unmapped",
            BamFlag::Reverse => "reverse",
            BamFlag::Read1 => "read-1",
            BamFlag::Read2 => "read-2",
            BamFlag::Secondary => "secondary",
            BamFlag::QCFail => "qc-fail",
            BamFlag::Duplicate => "duplicate",
            BamFlag::Supplementary => "supplementary",
        }
    }
}

#[cfg(feature = "cli")]
impl std::fmt::Display for BamFlag {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> Result<(), std::fmt::Error> {
        write!(f, "{}", self.as_str())
    }
}

fn add_to_flag(flags: &[BamFlag], flag: &mut u16) -> Result<(), Error> {
    for f in flags.iter() {
        *flag |= *f as u16;
    }

    Ok(())
}

impl ReadFilter {
    pub fn new(count_orphans: bool, excl_flags: &[BamFlag], incl_flags: &[BamFlag]) -> Result<Self, Error> {
        let mut s = Self {
            inc_flag: 0,
            exc_flag: 0,
            count_orphans,
        };

        s.add_incl_flags(incl_flags)?;
        s.add_excl_flags(excl_flags)?;
        Ok(s)
    }
    pub fn add_excl_flags(&mut self, flags: &[BamFlag]) -> Result<(), Error> {
        let mut flag: u16 = 0;
        add_to_flag(flags, &mut flag)?;
        self.exc_flag = flag;

        Ok(())
    }

    pub fn add_incl_flags(&mut self, flags: &[BamFlag]) -> Result<(), Error> {
        let mut flag: u16 = 0;
        add_to_flag(flags, &mut flag)?;
        self.inc_flag = flag;

        Ok(())
    }

    #[inline(always)]
    pub fn check_read(&mut self, read: &Record) -> bool {
        // let mut pass;

        // check if orphan pair
        if !self.count_orphans && (read.is_paired() && !read.is_proper_pair()) {
            return false;
        }

        if self.inc_flag > 0 && (self.inc_flag & read.inner.core.flag) == 0 {
            return false;
        }

        if self.exc_flag > 0 && (self.exc_flag & read.inner.core.flag) > 0 {
            return false;
        }

        true
    }
}

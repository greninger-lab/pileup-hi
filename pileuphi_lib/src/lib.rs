pub(crate) mod alignment;
pub(crate) mod bamio;
pub(crate) mod baq;
pub(crate) mod basedepth_string;
pub(crate) mod cigar_resolve;
pub(crate) mod engine;
pub(crate) mod errors;
pub(crate) mod jobqueue;
pub(crate) mod output;
pub(crate) mod overlap;
pub(crate) mod params;
pub(crate) mod pileup_iterator;
pub(crate) mod pileup_string;
pub(crate) mod position_queue;
pub(crate) mod read_buf;
pub(crate) mod read_filter;
pub(crate) mod refseq;
pub(crate) mod threading;
pub(crate) mod utils;

pub use engine::{PileupEngine, PileupSink, PileupStream};
pub use jobqueue::setup_exit_handler;
pub use output::PileupCoordinate;

pub mod error {
    pub use crate::errors::{Error, ErrorKind};
}

pub mod outputs {
    pub use crate::basedepth_string::BaseDepthString;
    pub use crate::pileup_string::PileupString;
}

pub mod param {
    pub use crate::params::{InputParams, PileupParams, STDOUT_ARG_STR};
    pub use crate::read_filter::BamFlag;
}

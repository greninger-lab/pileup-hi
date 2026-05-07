use crossbeam::channel::SendError;
use log::error;
use rust_htslib::errors::Error as HtslibError;

#[derive(Debug)]
/// Internal error type used by pileup-hi
pub struct Error {
    kind: ErrorKind,
}

impl Error {
    /// Handler for non-typed errors like anyhow used by dependencies. Instead of copying entire
    /// source string we log it once.
    pub fn from_generic(e: Box<dyn std::error::Error>, kind: ErrorKind) -> Self {
        error!("{e}");
        Self { kind }
    }

    pub fn kind(&self) -> &ErrorKind {
        &self.kind
    }
}

impl From<std::io::Error> for Error {
    fn from(e: std::io::Error) -> Self {
        Self {
            kind: ErrorKind::IOError(e),
        }
    }
}

impl From<std::num::TryFromIntError> for Error {
    fn from(e: std::num::TryFromIntError) -> Self {
        Self {
            kind: ErrorKind::AnomalousData(format!("{e}")),
        }
    }
}

impl From<std::num::ParseIntError> for Error {
    fn from(e: std::num::ParseIntError) -> Self {
        Self {
            kind: ErrorKind::AnomalousData(format!("{e}")),
        }
    }
}

impl<T> From<SendError<T>> for Error {
    fn from(e: SendError<T>) -> Self {
        error!("{e}");
        Self {
            kind: ErrorKind::ThreadingError("Threading failure"),
        }
    }
}

impl From<HtslibError> for Error {
    fn from(e: HtslibError) -> Self {
        Self {
            kind: ErrorKind::HTSLibError(e),
        }
    }
}

impl From<ErrorKind> for Error {
    fn from(kind: ErrorKind) -> Self {
        Self { kind }
    }
}

impl From<std::str::Utf8Error> for Error {
    fn from(_: std::str::Utf8Error) -> Self {
        Self {
            kind: ErrorKind::IOError(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Byte sequence not valid UTF-8",
            )),
        }
    }
}

impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.kind {
            ErrorKind::IOError(s) => {
                write!(f, "Bad IO operation: {s}")
            }

            ErrorKind::ThreadingError(s) => {
                write!(f, "Threading error: {s}")
            }

            ErrorKind::BamNotSortedByCoordinate(posa, posb) => {
                write!(
                    f,
                    "BAM not coordinate-sorted: read with pos {posa} comes before read with pos {posb}."
                )
            }

            ErrorKind::BamNotSortedByReference(tida, tidb) => {
                write!(
                    f,
                    "BAM not reference-sorted: read with reference index {tida} comes before read with pos {tidb}."
                )
            }

            ErrorKind::HTSLibError(e) => {
                write!(f, "rust_htslib internal error: {e}")
            }

            ErrorKind::RefSeqError(s) => {
                write!(f, "Error occurred when loading/retrieving refseq: {s}")
            }

            ErrorKind::MateOverlapFailed(name) => {
                write!(f, "Mate overlap correction failed with read {name}")
            }

            ErrorKind::BAQFailed(name) => {
                write!(f, "BAQ realignment failed with read {name}")
            }

            ErrorKind::AnomalousData(s) => {
                write!(f, "Anomalous (unusable) data detected: {s}")
            }

            ErrorKind::BadInputRegions(s) => {
                write!(f, "Couldn't parse genomic region list: {s}")
            }

            ErrorKind::UnknownBamFlag(s) => {
                write!(f, "Bad BAM flag: {s}")
            }
        }
    }
}

impl std::error::Error for Error {}

pub type ReadQNAME = String;

#[derive(Debug)]
pub enum ErrorKind {
    IOError(std::io::Error),
    ThreadingError(&'static str),
    HTSLibError(HtslibError),
    AnomalousData(String),
    BadInputRegions(String),
    RefSeqError(&'static str),
    BamNotSortedByCoordinate(i64, i64),
    BamNotSortedByReference(i32, i32),
    MateOverlapFailed(ReadQNAME),
    BAQFailed(ReadQNAME),
    UnknownBamFlag(Box<str>),
}

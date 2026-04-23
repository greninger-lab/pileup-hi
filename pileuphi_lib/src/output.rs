use crate::alignment::PileupAlignment;
use crate::errors::Error;
use crate::refseq::RefSeqHandle;
use crate::utils::OutputWriter;

/// The interface requirements for a pileup output. It needs to give ref information,
/// intake pileup alignments, update current ref info, display depth, and write itself.
pub trait OrderedPileupOutput: Send + Sync + Clone + std::fmt::Debug {
    /// Get the reference of the pileup
    #[allow(dead_code)]
    fn tid(&self) -> i32;

    /// Get the coordinate of the pileup
    #[allow(dead_code)]
    fn pos(&self) -> i64;

    /// Update internal data with pileup alignment
    fn intake(&mut self, p: &PileupAlignment, refseq: &RefSeqHandle) -> Result<(), Error>;
    /// Update reference data given ref num, pos, name, and sequence
    fn set_ref_info(&mut self, tid: i32, pos: i64, ref_name: &str, refseq: &RefSeqHandle);

    fn write<W: std::io::Write>(&mut self, writer: &mut W) -> Result<(), Error>;

    fn depth(&self) -> u32;

    fn clear(&mut self);

    #[allow(dead_code)]
    fn new() -> Self;
}

pub struct OutputFormat<T: OrderedPileupOutput> {
    output: T,
    writer: OutputWriter,
}

impl<T: OrderedPileupOutput> OutputFormat<T> {
    pub fn new(output: T, writer: OutputWriter) -> Self {
        Self { output, writer }
    }

    pub fn reject(&mut self) -> bool {
        self.output.clear();
        false
    }

    pub fn cur(&mut self) -> &mut T {
        &mut self.output
    }

    pub fn take(&mut self) -> Result<bool, Error> {
        self.output.write(&mut self.writer)?;
        Ok(true)
    }

    pub fn check(&mut self, emit: bool) -> Result<bool, Error> {
        if emit {
            self.take()
        } else {
            Ok(self.reject())
        }
    }
}

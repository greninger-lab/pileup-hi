use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    errors::{Error, ErrorKind},
    jobqueue::IntervalJobs,
    output::{OrderedPileupOutput, OutputFormat},
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::{RefSeq, RefSeqHandle},
    threading::ThreadPool,
    utils::get_writer_multi,
};

use log::{info, warn};
use std::sync::Arc;

pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
    dest: OutputDataDest,
    refseq: Option<RefSeq>,
}

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    pub fn initialize(in_params: InputParams, plp_params: PileupParams, output: T) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;
        let dest = OutputDataDest::from_string(&plp_params.output);

        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = in_params.region {
            create_region_queue(&region, header)?
        } else {
            intervals_from_header(header)?
        };

        let refseq = if let Some(ref fasta) = plp_params.refseq {
            if !std::fs::exists(std::path::Path::new(fasta))? {
                return Err(Error::from(ErrorKind::IOError(std::io::Error::new(
                    std::io::ErrorKind::NotFound,
                    "Fasta file {fasta} doesn't exist!",
                ))));
            }
            Some(RefSeq::new(fasta.clone()))
        } else {
            None
        };

        Ok(Self {
            intervals,
            plp_params,
            src,
            output,
            dest,
            refseq,
        })
    }

    pub fn run(self) -> Result<(), Error> {
        if self.intervals.is_empty() {
            return Ok(());
        }

        // remove old output file if it exists.
        if let OutputDataDest::File(ref f) = self.dest {
            if std::fs::exists(f)? {
                warn!("Output file {} already exists! Overwriting...", f);

                if let Err(e) = std::fs::remove_file(f) {
                    warn!("Failed to remove file {f}; {e}. Output will be appended...");
                };
            }
        }

        if self.src.has_index()? {
            info!("Found index for for input file {}", self.src.fname()?);
        }

        if self.plp_params.threads == 1 {
            self.run_single()
        } else if !self.src.has_index()? {
            warn!(
                "User asked for more than {} threads but file is unindexed. Running in single-thread mode...",
                self.plp_params.threads
            );
            self.run_single()
        } else {
            info!("Running with {} threads...", self.plp_params.threads);
            self.run_multi()
        }
    }

    fn get_refseq(&self, ref_name: &str) -> Result<RefSeqHandle, Error> {
        if let Some(ref refseq) = self.refseq {
            refseq.yield_handle(ref_name)
        } else {
            Ok(Arc::new(None))
        }
    }

    /// Use a single thread for both processing and writing.
    pub fn run_single(self) -> Result<(), Error> {
        for interval in self.intervals.iter() {
            let main_writer = get_writer_multi(&self.dest, BUFWRITER_CAP, true, false)?;

            let refseq_handle = self.get_refseq(&interval.name)?;

            let mut iterator = PileupIterator::new(
                &self.src,
                refseq_handle,
                &self.plp_params,
                OutputFormat::new(self.output.clone(), main_writer),
            )?;

            iterator.auto_loop2(interval)?;
        }
        Ok(())
    }

    /// Split up a list of input genomic intervals into smaller chunks to be processed in parallel. Chunks are first written to temporary output files before being merged into the user-specified output file.
    pub fn run_multi(self) -> Result<(), Error> {
        let mut jobs = IntervalJobs::new(
            &self.intervals,
            self.plp_params.coords_per_thread,
            self.plp_params.threads as i64,
            self.dest.clone(),
        );

        let mut pool = ThreadPool::new(self.plp_params.threads);
        let mut n_jobs = 0;

        while !jobs.is_completed() {
            jobs.merge_completed()?;

            if let Some(worker) = pool.get_available() {
                if let Some(job) = jobs.queue.pop_front() {
                    n_jobs += 1;

                    let refseq_handle = self.get_refseq(&job.interval.name)?;

                    worker.run(
                        n_jobs,
                        self.plp_params.clone(),
                        job,
                        self.src.clone(),
                        self.output.clone(),
                        refseq_handle,
                    );
                }
            }
        }
        jobs.conclude()
    }
}

use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    output::{
        generate_subfile_dests, OrderedPileupOutput, OutputFileMerge, OutputMethod, PileupOutputArray,
        FILE_MERGE_SINGLETON,
    },
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    refseq::RefSeq,
    utils::OutputWriter,
};

use std::sync::Arc;
use std::time::Instant;

const OUTPUT_ARRAY_YIELD_SIZE: i64 = 2000;
pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_BAM_READ_THREADS: usize = 2;

/// The default minimum number of coordinates to give each thread for processing.
/// This basically exists to prevent doing unnecessary work for very small regions.
/// Can be overridden if you need more horsepower for, say, high-depth regions.
pub const MIN_COORDS_PER_THREAD: i64 = 250_000;

use anyhow::Error;
use log::{info, warn};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::io::BufWriter;

pub struct PileupWorker {
    interval: GenomeInterval,
    params: PileupParams,
    src: BamDataSource,
}

impl PileupWorker {
    pub fn new(params: PileupParams, interval: GenomeInterval, src: BamDataSource) -> Self {
        Self { interval, params, src }
    }

    pub fn run<T>(&mut self, o: T, out: OutputWriter, refseq: RefSeqHandle)
    where
        T: OrderedPileupOutput + 'static,
    {
        let mut iterator = PileupIterator::new(
            &self.src,
            refseq,
            &self.params,
            o,
            OutputMethod::QueueForOutput(PileupOutputArray::new(
                std::cmp::min((self.interval.len() / 10).max(1), OUTPUT_ARRAY_YIELD_SIZE) as usize,
                out,
            )),
        )
        .unwrap();

        iterator.auto_loop2(&self.interval).unwrap();
    }
}

pub type RefSeqHandle = Option<Arc<Vec<u8>>>;

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
    dest: OutputDataDest,
    refseq: RefSeq,
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

        let refseq = if let Some(ref_file) = &plp_params.refseq {
            RefSeq::from_file(ref_file)?
        } else {
            RefSeq::blank()
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

    /// Use a single thread for both processing and writing.
    pub fn run_single(mut self) -> Result<(), Error> {
        for interval in self.intervals.iter() {
            self.refseq.load_seq(&interval.name)?;
            let lock = Box::new(BufWriter::with_capacity(2 * 1024 * 1024, std::io::stdout().lock()));

            let mut iterator = PileupIterator::new(
                &self.src,
                self.refseq.yield_handle(),
                &self.plp_params,
                self.output.clone(),
                OutputMethod::WriteDirectly(self.output.clone(), lock),
            )?;

            iterator.auto_loop2(interval)?;
        }
        Ok(())
    }

    /// Use separate threads for processing and writing. Each processing thread owns its IO readers for input BAM, index, and any other files.
    pub fn run_multi(mut self) -> Result<(), Error> {
        let outprefix = self.src.fname()?;

        let threadpool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.plp_params.threads)
            .build()
            .unwrap();

        for interval in &self.intervals {
            self.refseq.load_seq(&interval.name)?;
            let mut output_merge_lock = FILE_MERGE_SINGLETON.lock().expect("Failed to lock output file mutex");

            // we update the singleton tracking temp output files in case the program exits before finishing. This way the files
            // are marked for deletion during exit handling.
            *output_merge_lock = OutputFileMerge {
                outfile: self.dest.clone(),
                subfiles: generate_subfile_dests(&outprefix, self.plp_params.threads - 1, "temp.txt"),
            };

            // we use thread-local copy so we can drop the mutex lock
            let local_outputs = output_merge_lock.clone();
            drop(output_merge_lock);

            let per_thread_intervals =
                if self.plp_params.threads as i64 * self.plp_params.coords_per_thread > interval.len() {
                    interval.chunks(self.plp_params.coords_per_thread)
                } else {
                    interval.n_chunks(self.plp_params.threads as i64)
                }
                .collect::<Vec<GenomeInterval>>();

            info!(
                "Split ref {} into {} chunks...",
                interval.name,
                per_thread_intervals.len(),
            );

            let src = &self.src.clone();

            let before = Instant::now();

            let refhandle = self.refseq.yield_handle();

            threadpool.install(|| {
                per_thread_intervals.par_iter().enumerate().for_each(|(i, chunk)| {
                    let mut worker = PileupWorker::new(self.plp_params.clone(), chunk.clone(), src.clone());
                    let writer = local_outputs.get_writer(i).expect("failed to get writer");
                    worker.run(self.output.clone(), writer, refhandle.clone());
                });
            });

            let main_writer = local_outputs.get_writer(0)?;
            local_outputs.merge(main_writer)?;

            info!(
                "Tid {} completed in {} seconds...",
                interval.tid,
                before.elapsed().as_secs()
            );
        }

        Ok(())
    }
}

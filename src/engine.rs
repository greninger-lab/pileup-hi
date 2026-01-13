use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    output::{
        generate_subfile_dests, OrderedPileupOutput, OutputFileMerge, OutputMethod, PileupOutputArray,
        FILE_MERGE_SINGLETON,
    },
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    utils::{determine_thread_scheme, OutputWriter},
};

use std::time::Instant;

const OUTPUT_ARRAY_YIELD_SIZE: i64 = 2000;
pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_COORDS_PER_THREAD: usize = 250_000; // heuristic from benchmarking
pub const MIN_BAM_READ_THREADS: usize = 2;

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

    pub fn run<T>(&mut self, o: T, out: OutputWriter, read_threads: usize)
    where
        T: OrderedPileupOutput + 'static,
    {
        let mut iterator = PileupIterator::new(
            &self.src,
            std::slice::from_ref(&self.interval),
            &self.params,
            o,
            OutputMethod::QueueForOutput(PileupOutputArray::new(
                std::cmp::min((self.interval.len() / 10).max(1), OUTPUT_ARRAY_YIELD_SIZE) as usize,
                out,
            )),
            read_threads,
        )
        .unwrap();

        iterator.auto_loop().unwrap();
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
    dest: OutputDataDest,
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

        Ok(Self {
            intervals,
            plp_params,
            src,
            output,
            dest,
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
    pub fn run_single(self) -> Result<(), Error> {
        let lock = Box::new(BufWriter::with_capacity(2 * 1024 * 1024, std::io::stdout().lock()));
        let mut iterator = PileupIterator::new(
            &self.src,
            &self.intervals,
            &self.plp_params,
            self.output.clone(),
            OutputMethod::WriteDirectly(lock),
            MIN_BAM_READ_THREADS,
        )?;

        iterator.auto_loop()
        // iterator.auto_loop(&self.intervals[0], true)
    }

    /// Use separate threads for processing and writing. Each processing thread owns its IO readers for input BAM, index, and any other files.
    pub fn run_multi(self) -> Result<(), Error> {
        let outprefix = self.src.fname()?;

        let threadpool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.plp_params.threads)
            .build()
            .unwrap();

        for interval in &self.intervals {
            let thread_scheme = determine_thread_scheme(self.plp_params.threads, interval.len());

            let mut output_merge_lock = FILE_MERGE_SINGLETON.lock().expect("Failed to lock output file mutex");

            // we update the singleton tracking temp output files in case the program exits before finishing. This way the files
            // are marked for deletion during exit handling.
            *output_merge_lock = OutputFileMerge {
                outfile: self.dest.clone(),
                subfiles: generate_subfile_dests(&outprefix, thread_scheme.worker_threads - 1, "temp.txt"),
            };

            // we use thread-local copy so we can drop the mutex lock
            let local_outputs = output_merge_lock.clone();
            drop(output_merge_lock);

            let per_thread_intervals = interval
                .n_chunks(thread_scheme.worker_threads as i64)
                .collect::<Vec<GenomeInterval>>();

            info!(
                "Split ref {} into {} chunks...",
                interval.name,
                per_thread_intervals.len(),
            );

            let src = &self.src.clone();

            let before = Instant::now();

            threadpool.install(|| {
                per_thread_intervals.par_iter().enumerate().for_each(|(i, chunk)| {
                    let mut worker = PileupWorker::new(self.plp_params.clone(), chunk.clone(), src.clone());
                    let writer = local_outputs.get_writer(i).expect("failed to get writer");
                    worker.run(self.output.clone(), writer, thread_scheme.read_threads);
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

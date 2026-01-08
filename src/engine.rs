use crate::{
    bamio::{BamDataSource, BamReader, OutputDataDest},
    output::{
        generate_subfile_dests, OrderedPileupOutput, OutputFileMerge, OutputMethod, PileupOutputArray,
        TempOutputHandle, FILE_MERGE_SINGLETON,
    },
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
};

use std::time::Instant;

const OUTPUT_ARRAY_YIELD_SIZE: usize = 2000;
pub const BUFWRITER_CAP: usize = 2 * 1024 * 1024;
pub const MIN_COORDS_PER_THREAD: usize = 1000;

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

    pub fn run<T>(&mut self, o: T, out: TempOutputHandle, id: usize)
    where
        T: OrderedPileupOutput + 'static,
    {
        let mut iterator = PileupIterator::new(
            &self.src,
            std::slice::from_ref(&self.interval),
            &self.params,
            o,
            OutputMethod::QueueForOutput(
                PileupOutputArray::new(
                    1_000_000,
                    std::cmp::min((self.interval.len() / 10).max(1), OUTPUT_ARRAY_YIELD_SIZE),
                    id,
                    out,
                )
                .unwrap(),
            ),
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
                std::fs::remove_file(f)?;
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
            let mut output_merge_lock = FILE_MERGE_SINGLETON.lock().expect("Failed to lock output file mutex");

            let n_chunks = if interval.len() < MIN_COORDS_PER_THREAD {
                1
            } else {
                self.plp_params.threads as i64
            };

            // we update the singleton tracking temp output files in case the program exits before finishing. This way the files
            // are marked for deletion during exit handling.
            *output_merge_lock = OutputFileMerge {
                outfile: self.dest.clone(),
                subfiles: generate_subfile_dests(&outprefix, n_chunks as usize - 1, "temp.txt"),
            };

            // we use thread-local copy so we can drop the mutex lock
            let local_outputs = output_merge_lock.clone();
            drop(output_merge_lock);

            let per_thread_intervals = interval.n_chunks(n_chunks).collect::<Vec<GenomeInterval>>();

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
                    worker.run(self.output.clone(), writer, i);
                });
            });

            let main_writer = local_outputs.get_writer(0)?;
            local_outputs.merge(main_writer.writer)?;

            info!(
                "Tid {} completed in {} seconds...",
                interval.tid,
                before.elapsed().as_secs()
            );
        }

        Ok(())
    }
}

use crate::{
    bamio::{BamDataSource, BamReader},
    output::{
        merge_temp_outputs, OrderedPileupOutput, OutputMethod, PileupOutputArray, TempOutputHandle,
        TEMP_FILES,
    },
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    position_queue::{create_region_queue, intervals_from_header, GenomeInterval},
    utils::temp_fname,
};

const OUTPUT_ARRAY_YIELD_SIZE: usize = 2000;
const BUFWRITER_CAP: usize = 2 * 1024 * 1024;

use anyhow::Error;
use log::{info, warn};
use rayon::iter::{IndexedParallelIterator, IntoParallelRefIterator, ParallelIterator};
use std::io::{stdout, BufWriter};

pub struct PileupWorker {
    interval: GenomeInterval,
    params: PileupParams,
    src: BamDataSource,
}

pub struct DummyOutputWriter {}
impl std::io::Write for DummyOutputWriter {
    fn write(&mut self, _buf: &[u8]) -> std::io::Result<usize> {
        Ok(0)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(())
    }
}

impl PileupWorker {
    pub fn new(params: PileupParams, interval: GenomeInterval, src: BamDataSource) -> Self {
        Self {
            interval,
            params,
            src,
        }
    }

    pub fn run<T>(&mut self, o: T, out: TempOutputHandle, id: usize)
    where
        T: OrderedPileupOutput + 'static,
    {
        let iterator = PileupIterator::new(
            &self.src,
            &self.params,
            o,
            OutputMethod::<DummyOutputWriter, T>::QueueForOutput(
                PileupOutputArray::new(
                    1_000_000,
                    std::cmp::min(self.interval.len() / 10, OUTPUT_ARRAY_YIELD_SIZE),
                    id,
                    out,
                )
                .unwrap(),
            ),
        )
        .unwrap();

        iterator
            ._auto_loop_yield_batch(std::slice::from_ref(&self.interval))
            .unwrap()
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: Vec<GenomeInterval>,
    plp_params: PileupParams,
    src: BamDataSource,
    output: T,
}

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    pub fn initialize(
        in_params: InputParams,
        plp_params: PileupParams,
        output: T,
    ) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;

        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = &in_params.region {
            create_region_queue(region, header)?
        } else {
            intervals_from_header(header)?
        };

        Ok(Self {
            intervals,
            plp_params,
            src,
            output,
        })
    }

    pub fn run(self) -> Result<(), Error> {
        if self.plp_params.threads == 1 {
            self.run_single()
        } else if !self.src.has_index()? {
            warn!("User asked for more than {} threads but file is unindexed. Running in single-thread mode...", self.plp_params.threads);
            self.run_single()
        } else {
            info!("Running with {} threads...", self.plp_params.threads);
            self.run_multi()
        }
    }

    /// Use a single thread for both processing and writing.
    pub fn run_single(self) -> Result<(), Error> {
        for interval in self.intervals {
            let lock = BufWriter::with_capacity(2 * 1024 * 1024, std::io::stdout().lock());
            let mut iterator = PileupIterator::new(
                &self.src,
                &self.plp_params,
                self.output.clone(),
                OutputMethod::WriteDirectly(lock),
            )?;

            iterator._auto_loop_output_each(&[interval])?;
        }
        Ok(())
    }

    /// Use separate threads for processing and writing. Each processing thread owns its IO readers for input BAM, index, and any other files.
    pub fn run_multi(self) -> Result<(), Error> {
        let outprefix = self.src.fname()?;

        for interval in &self.intervals {
            let mut merge_map = Vec::with_capacity(self.plp_params.threads);
            merge_map.push("STDOUT".to_string());
            for i in 1..self.plp_params.threads {
                merge_map.push(temp_fname(&outprefix, &i.to_string(), "temp.txt"))
            }

            // we update TEMP_FILES in case the program exits before finishing. This way the files
            // are marked for deletion during exit handling.
            *TEMP_FILES.lock().expect("Failed to lock output file mutex") = merge_map.clone();

            let per_thread_intervals = interval
                .n_chunks(self.plp_params.threads as i64)
                .collect::<Vec<GenomeInterval>>();

            let threadpool = rayon::ThreadPoolBuilder::new()
                .num_threads(self.plp_params.threads)
                .build()
                .unwrap();

            info!(
                "Split tid {} into {} chunks for {} threads...",
                interval.tid,
                per_thread_intervals.len(),
                self.plp_params.threads
            );

            let src = &self.src.clone();

            threadpool.install(|| {
                per_thread_intervals
                    .par_iter()
                    .enumerate()
                    .for_each(|(i, chunk)| {
                        let mut worker =
                            PileupWorker::new(self.plp_params.clone(), chunk.clone(), src.clone());
                        worker.run(
                            self.output.clone(),
                            TempOutputHandle::new(&merge_map[i], BUFWRITER_CAP).unwrap(),
                            i,
                        );
                    });
            });

            info!(
                "Processing for tid {} completed. Deleting intermediate files...",
                interval.tid
            );

            merge_temp_outputs(
                &merge_map,
                BufWriter::with_capacity(BUFWRITER_CAP, stdout().lock()),
            )?;
        }

        Ok(())
    }
}

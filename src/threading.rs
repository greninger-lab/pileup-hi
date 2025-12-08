#![allow(dead_code, unused_imports)]
use crate::{
    bamio::{BamDataSource, BamReader},
    basedepth_string::BaseDepthString,
    output::{OrderedPileupOutput, OutputMethod, PileupOutputAggregator},
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    pileup_string::PileupString,
    position_queue::{create_region_queue, GenomeInterval, PositionQueue},
};

use rayon::prelude::*;
use rust_htslib::bam::pileup::Pileup;

const DEFAULT_READ_LEN: usize = 150;

use std::{collections::VecDeque, io::BufWriter};
use std::{ptr::with_exposed_provenance, thread::JoinHandle};

use anyhow::Error;
use crossbeam::channel::Sender;

pub enum PileupWorkerState {
    Off,
    Running(JoinHandle<()>),
}

pub struct PileupWorker {
    interval: GenomeInterval,
    id: usize,
    state: PileupWorkerState,
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
    pub fn new(p: PileupParams, interval: GenomeInterval, id: usize, src: BamDataSource) -> Self {
        Self {
            interval: interval.clone(),
            id,
            params: p.clone(),
            state: PileupWorkerState::Off,
            src: src.clone(),
        }
    }

    pub fn run<T>(&mut self, o: T, queue_handle: Sender<T>) -> Vec<T>
    where
        T: OrderedPileupOutput + 'static,
    {
        let iterator = PileupIterator::new(
            &self.src,
            &self.params,
            o,
            OutputMethod::<DummyOutputWriter, T>::QueueForOutput(queue_handle, Vec::with_capacity(10_000)),
        )
        .unwrap();

        iterator
            ._auto_loop_yield_batch(&PositionQueue {
                queue: vec![self.interval.clone()],
            })
            .unwrap()
    }

    pub fn wait(self) -> Result<(), Error> {
        match self.state {
            PileupWorkerState::Off => anyhow::bail!("Attempted to join a deactivated worker!"),
            PileupWorkerState::Running(j) => Ok(j.join().unwrap()),
        }
    }
}

pub struct PileupEngine<T: OrderedPileupOutput> {
    intervals: PositionQueue,
    read_size: usize,
    in_params: InputParams,
    plp_params: PileupParams,
    workers: Vec<PileupWorker>,
    src: BamDataSource,
    output: T,
}

impl<T: OrderedPileupOutput + 'static> PileupEngine<T> {
    pub fn initialize(in_params: InputParams, plp_params: PileupParams, output: T) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;
        let read_size = BamReader::sample_read_length(&src).unwrap_or(DEFAULT_READ_LEN);

        let tempreader = BamReader::new(&src, 1)?;
        let header = &tempreader.header;

        let intervals = if let Some(region) = &in_params.region {
            create_region_queue(region, header)?
        } else {
            PositionQueue::new(header)?
        };

        Ok(Self {
            intervals,
            workers: Vec::with_capacity(plp_params.threads),
            read_size,
            in_params,
            plp_params,
            src,
            output,
        })
    }

    pub fn run(self) -> Result<(), Error> {
        if self.intervals.len() == 1 || self.plp_params.threads == 1 {
            self.run_single()
        } else {
            self.run_multi()
        }
    }

    pub fn run_single(self) -> Result<(), Error> {
        let lock = std::io::stdout().lock();
        let mut iterator = PileupIterator::new(
            &self.src,
            &self.plp_params,
            self.output,
            OutputMethod::WriteDirectly(lock),
        )?;
        iterator._auto_loop(&self.intervals)
    }

    pub fn run_multi(self) -> Result<(), Error> {
        let threads = self.plp_params.threads;

        for interval in self.intervals.queue {
            let mut agg: PileupOutputAggregator<T> = PileupOutputAggregator::new();
            agg.run();
            let output_handle = agg.get_output_handle().unwrap();

            let subintervals = interval.chunks(1_000_000).collect::<VecDeque<GenomeInterval>>();

            eprintln!("Number of intervals: {}", subintervals.len());

            let threadpool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            let src = &self.src.clone();

            // thank you Seth Stadick for this this blazingly-fast rayon usage pattern.
            threadpool.install(|| {
                subintervals
                    .par_iter()
                    .enumerate()
                    .flat_map(|(i, chunk)| {
                        let mut worker = PileupWorker::new(self.plp_params.clone(), chunk.clone(), i, src.clone());
                        worker.run(self.output.clone(), output_handle.clone())
                    })
                    .for_each(|o| {
                        output_handle.send(o).unwrap();
                    });
            });
            drop(output_handle);
            agg.terminate()?;
        }

        Ok(())
    }
}

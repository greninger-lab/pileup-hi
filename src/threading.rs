#![allow(dead_code, unused_imports)]
use crate::{
    bamio::{BamDataSource, BamReader},
    output::PileupOutputAggregator,
    params::{InputParams, PileupParams},
    pileup_iterator::PileupIterator,
    pileup_writer::PileupString,
    position_queue::{create_region_queue, GenomeInterval, PositionQueue},
};

use std::collections::VecDeque;
use std::thread::JoinHandle;

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
    output_handle: Sender<PileupString>,
    src: BamDataSource,
}

impl PileupWorker {
    pub fn new(
        p: &PileupParams,
        interval: &GenomeInterval,
        id: usize,
        output_handle: Sender<PileupString>,
        src: &BamDataSource,
    ) -> Self {
        Self {
            interval: interval.clone(),
            id,
            params: p.clone(),
            state: PileupWorkerState::Off,
            output_handle,
            src: src.clone(),
        }
    }

    pub fn run(&mut self) {
        let p = self.params.clone();
        let i = self.interval.clone();
        let o = self.output_handle.clone();
        let s = self.src.clone();

        let j = std::thread::spawn(move || {
            let mut iterator = PileupIterator::new(&s, &p, Some(o)).unwrap();

            iterator
                ._auto_loop(&PositionQueue { queue: vec![i] })
                .unwrap();
        });

        self.state = PileupWorkerState::Running(j);
    }

    pub fn wait(self) -> Result<(), Error> {
        match self.state {
            PileupWorkerState::Off => anyhow::bail!("Attempted to join a deactivated worker!"),
            PileupWorkerState::Running(j) => Ok(j.join().unwrap()),
        }
    }
}

// A very simple driver of multiple concurrent pileup iterators.
pub struct PileupEngine {
    intervals: PositionQueue,
    read_size: usize,
    in_params: InputParams,
    plp_params: PileupParams,
    workers: Vec<PileupWorker>,
    src: BamDataSource,
}

impl PileupEngine {
    pub fn initialize(in_params: InputParams, plp_params: PileupParams) -> Result<Self, Error> {
        let src = BamDataSource::from_string(&in_params.file)?;
        let read_size = BamReader::sample_read_length(&src)?;

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
        })
    }

    pub fn run(&mut self) -> Result<(), Error> {
        self.run_single()
        // if self.intervals.len() > 1 {
        //     self.run_single()
        // } else {
        //     self.run_single()
        // }
    }

    pub fn run_single(&mut self) -> Result<(), Error> {
        let mut iterator = PileupIterator::new(&self.src, &self.plp_params, None)?;
        iterator._auto_loop(&self.intervals)
    }

    //     pub fn run_multi(&mut self) -> Result<(), Error> {
    //         for interval in &self.intervals.queue {
    //             let mut output: PileupOutputAggregator<PileupString> = PileupOutputAggregator::new();
    //             output.run();

    //             let mut subintervals: VecDeque<GenomeInterval> = interval
    //                 .n_chunks(self.plp_params.threads.try_into()?)
    //                 .collect();

    //             while !subintervals.is_empty() {
    //                 for i in 0..self.plp_params.threads {
    //                     if let Some(chunk) = subintervals.pop_front() {
    //                         self.workers.push(PileupWorker::new(
    //                             &self.plp_params,
    //                             &chunk,
    //                             i,
    //                             output.get_output_handle().unwrap(),
    //                         ));

    //                         self.workers[i].run();
    //                     }
    //                 }

    //                 for worker in self.workers.drain(..) {
    //                     worker.wait()?;
    //                 }
    //             }

    //             output.terminate()?;
    //         }

    //         Ok(())
    //     }
}

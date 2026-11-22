//! Streaming trajectory output during a run.
//!
//! [`FrameSinkObserver`] adapts the run loop's [`Observer`] hook to a
//! [`FrameSink`] from `tpt-dynamics-data`, so that frames are written as the
//! simulation runs (in the chunked, compressed format of that crate) rather than
//! collected in memory. The sink is closed when the run finishes. The simple
//! binary writer in the parent module ([`crate::output::TrajectoryWriter`]) is
//! unchanged.

use tpt_dynamics_data::trajectory::{Frame, FrameSink, TrajectoryError};

use crate::error::SimError;
use crate::output::{Observation, Observer};

/// Maps a trajectory-writer error onto the simulation error type.
fn map_err(e: TrajectoryError) -> SimError {
    match e {
        TrajectoryError::Io(io) => SimError::Io(io),
        other => SimError::Invalid(format!("trajectory sink: {other}")),
    }
}

/// An [`Observer`] that pushes one [`Frame`] per output step into a [`FrameSink`].
pub struct FrameSinkObserver<S: FrameSink + Send> {
    sink: S,
    n_atoms: usize,
    every: u64,
}

impl<S: FrameSink + Send> FrameSinkObserver<S> {
    /// Wraps `sink` for systems of `n_atoms` atoms, writing every `every` steps
    /// (values below 1 are treated as 1).
    pub fn new(sink: S, n_atoms: usize, every: u64) -> Self {
        Self { sink, n_atoms, every: every.max(1) }
    }

    /// Consumes the observer and returns the sink (after it has been closed by
    /// the run loop, or call [`Observer::finish`] first).
    pub fn into_sink(self) -> S {
        self.sink
    }
}

impl<S: FrameSink + Send> Observer for FrameSinkObserver<S> {
    fn interval(&self) -> u64 {
        self.every
    }

    fn observe(&mut self, obs: &Observation<'_>) -> Result<(), SimError> {
        let s = obs.system;
        if s.len() != self.n_atoms {
            return Err(SimError::Invalid(format!(
                "trajectory sink expects {} atoms, system has {}",
                self.n_atoms,
                s.len()
            )));
        }
        let l = s.sim_box.lengths;
        let frame = Frame {
            step: obs.energies.step,
            time: obs.energies.time,
            box_lengths: [l.x, l.y, l.z],
            positions: s.positions.iter().map(|p| [p.x, p.y, p.z]).collect(),
        };
        self.sink.push_frame(&frame).map_err(map_err)
    }

    fn finish(&mut self) -> Result<(), SimError> {
        self.sink.close().map_err(map_err)
    }
}

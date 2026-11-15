//! Streaming trajectory writer: frames written during a run through the
//! tpt-dynamics-data `FrameSink` are readable and match the positions the run
//! produced at the same steps.

use std::io::{self, Cursor, Write};
use std::sync::{Arc, Mutex};

use tpt_dynamics_core::Vec3;
use tpt_dynamics_data::trajectory::{Precision, TrajectoryReader, TrajectoryWriter};
use tpt_dynamics_sim::error::SimError;
use tpt_dynamics_sim::forcefield::{LennardJones, Modifier};
use tpt_dynamics_sim::integrator::{ForceEvaluator, VelocityVerlet};
use tpt_dynamics_sim::output::sink::FrameSinkObserver;
use tpt_dynamics_sim::output::{Observation, Observer};
use tpt_dynamics_sim::pbc::SimBox;
use tpt_dynamics_sim::simulation::{maxwell_boltzmann_velocities, Simulation};
use tpt_dynamics_sim::system::{NewAtom, System, UnitSystem};

/// A writer whose bytes remain reachable after the observer is dropped.
#[derive(Clone)]
struct SharedBuf(Arc<Mutex<Vec<u8>>>);

impl Write for SharedBuf {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        self.0.lock().unwrap().extend_from_slice(buf);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

type Seen = Arc<Mutex<Vec<(u64, Vec<[f64; 3]>)>>>;

/// Records the step and positions seen at each output step.
struct Capture {
    every: u64,
    seen: Seen,
}

impl Observer for Capture {
    fn interval(&self) -> u64 {
        self.every
    }

    fn observe(&mut self, o: &Observation<'_>) -> Result<(), SimError> {
        let pos = o.system.positions.iter().map(|p| [p.x, p.y, p.z]).collect();
        self.seen.lock().unwrap().push((o.energies.step, pos));
        Ok(())
    }
}

fn lj_cubic(n: usize, density: f64) -> System {
    let a = (1.0 / density).cbrt();
    let l = a * n as f64;
    let mut sys = System::new(SimBox::new(l, l, l).unwrap());
    for ix in 0..n {
        for iy in 0..n {
            for iz in 0..n {
                let p =
                    Vec3::new((ix as f64 + 0.5) * a, (iy as f64 + 0.5) * a, (iz as f64 + 0.5) * a);
                sys.push_atom(NewAtom {
                    name: "AR".into(),
                    residue_name: "LJ".into(),
                    residue_id: 1,
                    element: "Ar".into(),
                    position: p,
                    velocity: Vec3::ZERO,
                    mass: 1.0,
                    charge: 0.0,
                    atom_type: 0,
                });
            }
        }
    }
    sys
}

#[test]
fn frames_written_during_a_run_match_the_positions() {
    let mut sys = lj_cubic(6, 0.8);
    let n = sys.len();
    maxwell_boltzmann_velocities(&mut sys, 1.0, 3, UnitSystem::Reduced);
    let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
    let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
    let mut sim =
        Simulation::new(sys, eval, Box::new(VelocityVerlet), None, 0.002, UnitSystem::Reduced)
            .unwrap();

    let bytes = Arc::new(Mutex::new(Vec::new()));
    let writer = TrajectoryWriter::new(SharedBuf(bytes.clone()), n, Precision::F64, 2).unwrap();
    let seen: Seen = Arc::new(Mutex::new(Vec::new()));
    let mut observers: Vec<Box<dyn Observer>> = vec![
        Box::new(FrameSinkObserver::new(writer, n, 5)),
        Box::new(Capture { every: 5, seen: seen.clone() }),
    ];
    sim.run(23, &mut observers).unwrap();
    drop(observers);

    let captured = seen.lock().unwrap().clone();
    assert_eq!(captured.len(), 5, "steps 0, 5, 10, 15, 20");
    let data = bytes.lock().unwrap().clone();
    let mut reader = TrajectoryReader::open(Cursor::new(data)).unwrap();
    assert_eq!(reader.n_frames(), captured.len() as u64);
    assert_eq!(reader.n_atoms(), n);
    for (i, (step, pos)) in captured.iter().enumerate() {
        let frame = reader.frame(i as u64).unwrap();
        assert_eq!(frame.step, *step);
        assert_eq!(frame.positions.len(), n);
        for (a, b) in frame.positions.iter().zip(pos) {
            for k in 0..3 {
                assert_eq!(a[k].to_bits(), b[k].to_bits(), "frame {i}");
            }
        }
    }
}

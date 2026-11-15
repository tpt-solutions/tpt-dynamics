//! Run-time output: the [`Observer`] hook, a plain-binary trajectory format,
//! and a CSV energy log.
//!
//! # Trajectory layout (version 1, little-endian)
//!
//! | Field        | Type            | Notes                                   |
//! |--------------|-----------------|-----------------------------------------|
//! | `magic`        | 8 bytes         | ASCII `TPTTRJ01`                        |
//! | `version`      | `u32`           | `1`                                     |
//! | `n_atoms`    | `u64`           | atoms per frame, constant for the file  |
//! | frames       | see below       | frames repeat until end of file         |
//!
//! Each frame:
//!
//! | Field     | Type                | Notes                              |
//! |-----------|---------------------|------------------------------------|
//! | step      | `u64`               | integer step index                 |
//! | time      | `f64`               | simulation time in ps              |
//! | box       | 3 × `f64`           | orthorhombic edge lengths in nm    |
//! | positions | `n_atoms` × 3 × `f32` | x, y, z in nm, atom order as system |
//!
//! Positions are stored in single precision to halve the file size; this is
//! sufficient for analysis but not for restarting a run (checkpointing is a
//! later phase).

pub mod sink;

use std::io::{Read, Write};

use crate::error::SimError;
use crate::system::System;

/// Magic bytes at the start of a trajectory file.
pub const TRAJ_MAGIC: &[u8; 8] = b"TPTTRJ01";
/// Trajectory format version written by [`TrajectoryWriter`].
pub const TRAJ_VERSION: u32 = 1;

/// Scalar observables for one step.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Energies {
    /// Step index.
    pub step: u64,
    /// Simulation time in ps.
    pub time: f64,
    /// Kinetic energy in kJ/mol.
    pub kinetic: f64,
    /// Potential energy in kJ/mol.
    pub potential: f64,
    /// `kinetic + potential` in kJ/mol.
    pub total: f64,
    /// Instantaneous temperature in K.
    pub temperature: f64,
    /// Pressure in kJ mol⁻¹ nm⁻³ (multiply by 16.6054 for bar).
    pub pressure: f64,
}

/// What an observer sees at an output step.
#[derive(Clone, Copy, Debug)]
pub struct Observation<'a> {
    /// The system state (positions, velocities, box).
    pub system: &'a System,
    /// Scalar observables for the same state.
    pub energies: Energies,
}

/// A hook called by the run loop at each output step.
pub trait Observer: Send {
    /// Call interval in steps (values below 1 are treated as 1).
    fn interval(&self) -> u64;

    /// Receives one observation.
    ///
    /// # Errors
    /// Returns an error if the output cannot be written.
    fn observe(&mut self, observation: &Observation<'_>) -> Result<(), SimError>;

    /// Called once at the end of a run; flushes buffers.
    ///
    /// # Errors
    /// Returns an error if the flush fails.
    fn finish(&mut self) -> Result<(), SimError> {
        Ok(())
    }
}

/// Writes frames in the binary trajectory format described in the module docs.
pub struct TrajectoryWriter<W: Write> {
    writer: W,
    n_atoms: usize,
    every: u64,
}

impl<W: Write> TrajectoryWriter<W> {
    /// Writes the file header and returns a writer for `n_atoms` atoms,
    /// called every `every` steps.
    ///
    /// # Errors
    /// Returns [`SimError::Io`] if the header cannot be written.
    pub fn new(mut writer: W, n_atoms: usize, every: u64) -> Result<Self, SimError> {
        writer.write_all(TRAJ_MAGIC)?;
        writer.write_all(&TRAJ_VERSION.to_le_bytes())?;
        writer.write_all(&(n_atoms as u64).to_le_bytes())?;
        Ok(Self { writer, n_atoms, every })
    }
}

impl<W: Write + Send> Observer for TrajectoryWriter<W> {
    fn interval(&self) -> u64 {
        self.every.max(1)
    }

    fn observe(&mut self, obs: &Observation<'_>) -> Result<(), SimError> {
        let s = obs.system;
        if s.len() != self.n_atoms {
            return Err(SimError::Invalid(format!(
                "trajectory expects {} atoms, system has {}",
                self.n_atoms,
                s.len()
            )));
        }
        self.writer.write_all(&obs.energies.step.to_le_bytes())?;
        self.writer.write_all(&obs.energies.time.to_le_bytes())?;
        let l = s.sim_box.lengths;
        for v in [l.x, l.y, l.z] {
            self.writer.write_all(&v.to_le_bytes())?;
        }
        for p in &s.positions {
            for c in [p.x, p.y, p.z] {
                self.writer.write_all(&(c as f32).to_le_bytes())?;
            }
        }
        Ok(())
    }

    fn finish(&mut self) -> Result<(), SimError> {
        self.writer.flush()?;
        Ok(())
    }
}

/// One frame read back from a trajectory file.
#[derive(Clone, Debug, PartialEq)]
pub struct TrajectoryFrame {
    /// Step index.
    pub step: u64,
    /// Time in ps.
    pub time: f64,
    /// Box edge lengths in nm.
    pub box_lengths: [f64; 3],
    /// Positions in nm.
    pub positions: Vec<[f32; 3]>,
}

/// A complete trajectory read from a stream.
#[derive(Clone, Debug, PartialEq)]
pub struct Trajectory {
    /// Atoms per frame.
    pub n_atoms: usize,
    /// Frames in file order.
    pub frames: Vec<TrajectoryFrame>,
}

fn read_exact_or_eof<R: Read>(r: &mut R, buf: &mut [u8]) -> Result<bool, SimError> {
    let mut filled = 0;
    while filled < buf.len() {
        let n = r.read(&mut buf[filled..])?;
        if n == 0 {
            if filled == 0 {
                return Ok(false);
            }
            return Err(SimError::parse("TRJ", 0, "truncated frame"));
        }
        filled += n;
    }
    Ok(true)
}

/// Reads a whole trajectory into memory.
///
/// # Errors
/// Returns [`SimError::Parse`] for a bad magic, unsupported version or a
/// truncated frame, and [`SimError::Io`] on read failure.
pub fn read_trajectory<R: Read>(mut reader: R) -> Result<Trajectory, SimError> {
    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic).map_err(|_| SimError::parse("TRJ", 0, "missing header"))?;
    if &magic != TRAJ_MAGIC {
        return Err(SimError::parse("TRJ", 0, "bad magic"));
    }
    let mut u32buf = [0u8; 4];
    reader.read_exact(&mut u32buf).map_err(|_| SimError::parse("TRJ", 0, "missing header"))?;
    if u32::from_le_bytes(u32buf) != TRAJ_VERSION {
        return Err(SimError::parse("TRJ", 0, "unsupported version"));
    }
    let mut u64buf = [0u8; 8];
    reader.read_exact(&mut u64buf).map_err(|_| SimError::parse("TRJ", 0, "missing header"))?;
    let n_atoms = u64::from_le_bytes(u64buf) as usize;

    let mut frames = Vec::new();
    loop {
        let mut head = [0u8; 8];
        if !read_exact_or_eof(&mut reader, &mut head)? {
            break;
        }
        let step = u64::from_le_bytes(head);
        let mut f8 = [0u8; 8];
        let mut get_f64 = |r: &mut R| -> Result<f64, SimError> {
            if !read_exact_or_eof(r, &mut f8)? {
                return Err(SimError::parse("TRJ", 0, "truncated frame"));
            }
            Ok(f64::from_le_bytes(f8))
        };
        let time = get_f64(&mut reader)?;
        let box_lengths = [get_f64(&mut reader)?, get_f64(&mut reader)?, get_f64(&mut reader)?];
        let mut positions = Vec::with_capacity(n_atoms.min(1 << 24));
        let mut c = [0u8; 4];
        for _ in 0..n_atoms {
            let mut p = [0f32; 3];
            for slot in &mut p {
                if !read_exact_or_eof(&mut reader, &mut c)? {
                    return Err(SimError::parse("TRJ", 0, "truncated frame"));
                }
                *slot = f32::from_le_bytes(c);
            }
            positions.push(p);
        }
        frames.push(TrajectoryFrame { step, time, box_lengths, positions });
    }
    Ok(Trajectory { n_atoms, frames })
}

/// Writes one CSV row per output step: `step,time,kinetic,potential,total,temperature,pressure`.
pub struct EnergyLogger<W: Write> {
    writer: W,
    every: u64,
}

impl<W: Write> EnergyLogger<W> {
    /// Writes the header line and returns a logger called every `every` steps.
    ///
    /// # Errors
    /// Returns [`SimError::Io`] if the header cannot be written.
    pub fn new(mut writer: W, every: u64) -> Result<Self, SimError> {
        writeln!(writer, "step,time,kinetic,potential,total,temperature,pressure")?;
        Ok(Self { writer, every })
    }
}

impl<W: Write + Send> Observer for EnergyLogger<W> {
    fn interval(&self) -> u64 {
        self.every.max(1)
    }

    fn observe(&mut self, obs: &Observation<'_>) -> Result<(), SimError> {
        let e = &obs.energies;
        writeln!(
            self.writer,
            "{},{:.6},{:.8e},{:.8e},{:.8e},{:.6},{:.8e}",
            e.step, e.time, e.kinetic, e.potential, e.total, e.temperature, e.pressure
        )?;
        Ok(())
    }

    fn finish(&mut self) -> Result<(), SimError> {
        self.writer.flush()?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pbc::SimBox;
    use crate::system::NewAtom;
    use tpt_dynamics_core::Vec3;

    fn tiny_system(n: usize) -> System {
        let mut s = System::new(SimBox::new(3.0, 3.0, 3.0).unwrap());
        for i in 0..n {
            s.push_atom(NewAtom {
                name: "AR".into(),
                residue_name: "LJ".into(),
                residue_id: 1,
                element: "Ar".into(),
                position: Vec3::new(0.25 * i as f64, 0.5, 1.0 + 0.1 * i as f64),
                velocity: Vec3::ZERO,
                mass: 1.0,
                charge: 0.0,
                atom_type: 0,
            });
        }
        s
    }

    #[test]
    fn trajectory_round_trip() {
        let sys = tiny_system(3);
        let mut w = TrajectoryWriter::new(Vec::new(), 3, 1).unwrap();
        for step in 0..4u64 {
            let e = Energies { step, time: 0.01 * step as f64, ..Default::default() };
            w.observe(&Observation { system: &sys, energies: e }).unwrap();
        }
        w.finish().unwrap();
        let bytes = w.writer.clone();
        assert_eq!(&bytes[..8], TRAJ_MAGIC);
        assert_eq!(bytes.len(), 8 + 4 + 8 + 4 * (8 + 8 + 24 + 3 * 12));
        let traj = read_trajectory(&bytes[..]).unwrap();
        assert_eq!(traj.n_atoms, 3);
        assert_eq!(traj.frames.len(), 4);
        assert_eq!(traj.frames[2].step, 2);
        assert!((traj.frames[2].time - 0.02).abs() < 1e-15);
        assert_eq!(traj.frames[0].box_lengths, [3.0, 3.0, 3.0]);
        assert!((traj.frames[0].positions[1][0] - 0.25).abs() < 1e-6);
    }

    #[test]
    fn trajectory_reader_rejects_bad_input() {
        assert!(read_trajectory(&b""[..]).is_err());
        assert!(read_trajectory(&b"NOTTRAJ!\x01\0\0\0\0\0\0\0\0\0\0\0"[..]).is_err());
        let sys = tiny_system(2);
        let mut w = TrajectoryWriter::new(Vec::new(), 2, 1).unwrap();
        w.observe(&Observation { system: &sys, energies: Energies::default() }).unwrap();
        let bytes = w.writer;
        for n in 0..bytes.len() {
            let r = read_trajectory(&bytes[..n]);
            if n == 20 {
                assert!(r.is_ok(), "header only is a valid empty trajectory");
            }
        }
        assert!(read_trajectory(&bytes[..bytes.len() - 1]).is_err());
    }

    #[test]
    fn trajectory_writer_rejects_wrong_atom_count() {
        let sys = tiny_system(2);
        let mut w = TrajectoryWriter::new(Vec::new(), 5, 1).unwrap();
        assert!(w.observe(&Observation { system: &sys, energies: Energies::default() }).is_err());
    }

    #[test]
    fn energy_log_has_header_and_rows() {
        let sys = tiny_system(2);
        let mut log = EnergyLogger::new(Vec::new(), 10).unwrap();
        assert_eq!(log.interval(), 10);
        log.observe(&Observation {
            system: &sys,
            energies: Energies {
                step: 10,
                time: 0.1,
                kinetic: 1.5,
                potential: -2.0,
                total: -0.5,
                temperature: 1.0,
                pressure: 0.3,
            },
        })
        .unwrap();
        let text = String::from_utf8(log.writer).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], "step,time,kinetic,potential,total,temperature,pressure");
        assert!(lines[1].starts_with("10,0.100000,"));
        assert_eq!(lines[1].split(',').count(), 7);
    }
}

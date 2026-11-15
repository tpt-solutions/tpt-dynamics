//! Checkpoint and restart: a versioned, little-endian binary snapshot of the
//! full dynamical state, protected by a CRC-32 trailer.
//!
//! Layout (version 1):
//!
//! ```text
//! magic        8 bytes   b"TPTCKPT\0"
//! version      u32       1
//! n_atoms      u64
//! step         u64       completed steps
//! dt           f64       time step the run was using (checked on restore)
//! box          6 × f64   lengths lx, ly, lz and tilts xy, xz, yz
//! positions    3n × f64
//! velocities   3n × f64  (leap-frog: half-step velocities)
//! thermostat   u32 count, then count × u64 words
//! barostat     u32 count, then count × u64 words
//! crc32        u32       CRC-32 of every preceding byte
//! ```
//!
//! Forces, the neighbour list and the constraint geometry are not stored; they
//! are recomputed from positions on restore. [`Simulation::checkpoint_bytes`]
//! rebuilds the neighbour list before encoding, so that a restored run and an
//! uninterrupted run follow the same floating-point operations and produce
//! bitwise-identical states (see the tests in this module).
//!
//! The decoder checks every length against the remaining input before
//! allocating, so truncated or corrupted input yields an error rather than a
//! panic or an oversized allocation.

use tpt_dynamics_core::Vec3;

use crate::crc32::crc32;
use crate::error::SimError;
use crate::pbc::SimBox;
use crate::simulation::Simulation;

/// File magic.
pub const MAGIC: &[u8; 8] = b"TPTCKPT\0";
/// Format version written by this crate.
pub const VERSION: u32 = 1;

/// Largest number of state words accepted for a thermostat or barostat. Real
/// implementations use at most a few words.
const MAX_STATE_WORDS: usize = 64;

/// Decoded checkpoint contents.
#[derive(Clone, Debug, PartialEq)]
pub struct Checkpoint {
    /// Number of atoms.
    pub n_atoms: usize,
    /// Completed steps.
    pub step: u64,
    /// Time step (ps).
    pub dt: f64,
    /// Periodic box.
    pub sim_box: SimBox,
    /// Positions (nm).
    pub positions: Vec<Vec3>,
    /// Velocities (nm/ps).
    pub velocities: Vec<Vec3>,
    /// Thermostat state words.
    pub thermostat: Vec<u64>,
    /// Barostat state words.
    pub barostat: Vec<u64>,
}

/// Serialises the simulation state. Use [`Simulation::checkpoint_bytes`], which
/// also refreshes the neighbour list.
pub fn encode(sim: &Simulation) -> Vec<u8> {
    let n = sim.system.len();
    let mut out = Vec::with_capacity(96 + 48 * n + 64);
    out.extend_from_slice(MAGIC);
    out.extend_from_slice(&VERSION.to_le_bytes());
    out.extend_from_slice(&(n as u64).to_le_bytes());
    out.extend_from_slice(&sim.step.to_le_bytes());
    out.extend_from_slice(&sim.dt.to_le_bytes());
    let b = sim.system.sim_box;
    for v in [b.lengths.x, b.lengths.y, b.lengths.z, b.tilt[0], b.tilt[1], b.tilt[2]] {
        out.extend_from_slice(&v.to_le_bytes());
    }
    for p in sim.system.positions.iter().chain(&sim.system.velocities) {
        out.extend_from_slice(&p.x.to_le_bytes());
        out.extend_from_slice(&p.y.to_le_bytes());
        out.extend_from_slice(&p.z.to_le_bytes());
    }
    for words in [
        sim.thermostat.as_ref().map(|t| t.save_state()).unwrap_or_default(),
        sim.barostat.as_ref().map(|b| b.save_state()).unwrap_or_default(),
    ] {
        out.extend_from_slice(&(words.len() as u32).to_le_bytes());
        for w in words {
            out.extend_from_slice(&w.to_le_bytes());
        }
    }
    let crc = crc32(&out);
    out.extend_from_slice(&crc.to_le_bytes());
    out
}

/// Bounds-checked little-endian reader.
struct Reader<'a> {
    data: &'a [u8],
    pos: usize,
}

impl<'a> Reader<'a> {
    fn take(&mut self, n: usize) -> Result<&'a [u8], SimError> {
        let end = self.pos.checked_add(n).ok_or_else(truncated)?;
        if end > self.data.len() {
            return Err(truncated());
        }
        let s = &self.data[self.pos..end];
        self.pos = end;
        Ok(s)
    }

    fn u32(&mut self) -> Result<u32, SimError> {
        let b = self.take(4)?;
        Ok(u32::from_le_bytes([b[0], b[1], b[2], b[3]]))
    }

    fn u64(&mut self) -> Result<u64, SimError> {
        let b = self.take(8)?;
        let mut a = [0u8; 8];
        a.copy_from_slice(b);
        Ok(u64::from_le_bytes(a))
    }

    fn f64(&mut self) -> Result<f64, SimError> {
        Ok(f64::from_bits(self.u64()?))
    }

    fn remaining(&self) -> usize {
        self.data.len() - self.pos
    }
}

fn truncated() -> SimError {
    SimError::Invalid("checkpoint is truncated".into())
}

impl Checkpoint {
    /// Parses and validates a checkpoint byte stream.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] for a wrong magic, unsupported version, CRC
    /// mismatch, truncation, inconsistent lengths, non-finite values or an
    /// invalid box.
    pub fn decode(bytes: &[u8]) -> Result<Self, SimError> {
        if bytes.len() < MAGIC.len() + 4 + 4 {
            return Err(truncated());
        }
        let (body, trailer) = bytes.split_at(bytes.len() - 4);
        let stored = u32::from_le_bytes([trailer[0], trailer[1], trailer[2], trailer[3]]);
        if crc32(body) != stored {
            return Err(SimError::Invalid("checkpoint CRC mismatch".into()));
        }
        Self::decode_body(body)
    }

    /// Parses the body (everything before the CRC trailer). Public to the crate
    /// so that tests can exercise the parser past the checksum.
    pub(crate) fn decode_body(body: &[u8]) -> Result<Self, SimError> {
        let mut r = Reader { data: body, pos: 0 };
        if r.take(MAGIC.len())? != MAGIC {
            return Err(SimError::Invalid("not a tpt-dynamics checkpoint".into()));
        }
        let version = r.u32()?;
        if version != VERSION {
            return Err(SimError::Invalid(format!("unsupported checkpoint version {version}")));
        }
        let n_atoms = usize::try_from(r.u64()?)
            .map_err(|_| SimError::Invalid("atom count does not fit in memory".into()))?;
        let step = r.u64()?;
        let dt = r.f64()?;
        let lengths = [r.f64()?, r.f64()?, r.f64()?];
        let tilt = [r.f64()?, r.f64()?, r.f64()?];
        let sim_box =
            SimBox::triclinic(lengths[0], lengths[1], lengths[2], tilt[0], tilt[1], tilt[2])?;
        // Each atom needs 48 bytes (positions and velocities). Check before allocating.
        let needed = n_atoms
            .checked_mul(48)
            .ok_or_else(|| SimError::Invalid("atom count overflows".into()))?;
        if r.remaining() < needed {
            return Err(truncated());
        }
        let read_vec = |r: &mut Reader<'_>| -> Result<Vec<Vec3>, SimError> {
            let mut v = Vec::with_capacity(n_atoms);
            for _ in 0..n_atoms {
                v.push(Vec3::new(r.f64()?, r.f64()?, r.f64()?));
            }
            Ok(v)
        };
        let positions = read_vec(&mut r)?;
        let velocities = read_vec(&mut r)?;
        if positions
            .iter()
            .chain(&velocities)
            .any(|p| !(p.x.is_finite() && p.y.is_finite() && p.z.is_finite()))
        {
            return Err(SimError::Invalid("checkpoint contains non-finite values".into()));
        }
        if !(dt.is_finite() && dt > 0.0) {
            return Err(SimError::Invalid("checkpoint time step is invalid".into()));
        }
        let read_words = |r: &mut Reader<'_>| -> Result<Vec<u64>, SimError> {
            let count = r.u32()? as usize;
            if count > MAX_STATE_WORDS || r.remaining() < count * 8 {
                return Err(SimError::Invalid("checkpoint state block is invalid".into()));
            }
            (0..count).map(|_| r.u64()).collect()
        };
        let thermostat = read_words(&mut r)?;
        let barostat = read_words(&mut r)?;
        if r.remaining() != 0 {
            return Err(SimError::Invalid("checkpoint has trailing bytes".into()));
        }
        Ok(Self { n_atoms, step, dt, sim_box, positions, velocities, thermostat, barostat })
    }

    /// Writes this state into `sim`, which must match the checkpoint.
    ///
    /// # Errors
    /// Returns [`SimError::Invalid`] if the atom count, time step or thermostat /
    /// barostat state layout does not match the simulation.
    pub fn apply(&self, sim: &mut Simulation) -> Result<(), SimError> {
        if self.n_atoms != sim.system.len() {
            return Err(SimError::Invalid(format!(
                "checkpoint has {} atoms, simulation has {}",
                self.n_atoms,
                sim.system.len()
            )));
        }
        if self.dt.to_bits() != sim.dt.to_bits() {
            return Err(SimError::Invalid(format!(
                "checkpoint time step {} differs from simulation time step {}",
                self.dt, sim.dt
            )));
        }
        if let Some(t) = sim.thermostat.as_mut() {
            t.load_state(&self.thermostat)?;
        } else if !self.thermostat.is_empty() {
            return Err(SimError::Invalid(
                "checkpoint has thermostat state but none is configured".into(),
            ));
        }
        if let Some(b) = sim.barostat.as_mut() {
            b.load_state(&self.barostat)?;
        } else if !self.barostat.is_empty() {
            return Err(SimError::Invalid(
                "checkpoint has barostat state but none is configured".into(),
            ));
        }
        sim.system.positions.clone_from(&self.positions);
        sim.system.velocities.clone_from(&self.velocities);
        sim.system.sim_box = self.sim_box;
        sim.step = self.step;
        sim.eval.invalidate_neighbours();
        sim.eval.evaluate(&mut sim.system)?;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{LennardJones, Modifier};
    use crate::integrator::{ForceEvaluator, LeapFrog, VelocityVerlet};
    use crate::rng::Rng;
    use crate::simulation::maxwell_boltzmann_velocities;
    use crate::system::UnitSystem;
    use crate::thermostat::VRescale;

    /// Builds a deterministic LJ + v-rescale simulation.
    fn build(integrator_leap_frog: bool) -> Simulation {
        let mut sys = crate::integrator::tests::lj_gas(5, 0.7, 31);
        maxwell_boltzmann_velocities(&mut sys, 0.8, 4, UnitSystem::Reduced);
        let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap();
        let eval = ForceEvaluator::new(Box::new(ff), 0.3).unwrap();
        let integ: Box<dyn crate::integrator::Integrator> =
            if integrator_leap_frog { Box::new(LeapFrog) } else { Box::new(VelocityVerlet) };
        let th = VRescale::new(0.8, 0.2, UnitSystem::Reduced, 77).unwrap();
        Simulation::new(sys, eval, integ, Some(Box::new(th)), 0.002, UnitSystem::Reduced).unwrap()
    }

    fn snapshot(sim: &Simulation) -> (Vec<Vec3>, Vec<Vec3>, u64) {
        (sim.system.positions.clone(), sim.system.velocities.clone(), sim.step)
    }

    #[test]
    fn crc_protects_against_corruption() {
        let mut sim = build(false);
        sim.run(3, &mut []).unwrap();
        let bytes = sim.checkpoint_bytes().unwrap();
        let mut bad = bytes.clone();
        let mid = bad.len() / 2;
        bad[mid] ^= 0x10;
        assert!(Checkpoint::decode(&bad).is_err());
        assert!(Checkpoint::decode(&bytes).is_ok());
    }

    #[test]
    fn restart_reproduces_uninterrupted_run_bitwise() {
        for leap in [false, true] {
            // Uninterrupted: 40 steps, checkpoint, 40 more.
            let mut a = build(leap);
            a.run(40, &mut []).unwrap();
            let bytes = a.checkpoint_bytes().unwrap();
            a.run(40, &mut []).unwrap();
            let reference = snapshot(&a);

            // Restarted: fresh simulation, restore, 40 more.
            let mut b = build(leap);
            b.restore_checkpoint(&bytes).unwrap();
            assert_eq!(b.step, 40);
            b.run(40, &mut []).unwrap();
            let restarted = snapshot(&b);

            assert_eq!(reference.2, restarted.2);
            for (x, y) in reference.0.iter().zip(&restarted.0) {
                assert_eq!(x.x.to_bits(), y.x.to_bits());
                assert_eq!(x.y.to_bits(), y.y.to_bits());
                assert_eq!(x.z.to_bits(), y.z.to_bits());
            }
            for (x, y) in reference.1.iter().zip(&restarted.1) {
                assert_eq!(x.x.to_bits(), y.x.to_bits());
                assert_eq!(x.y.to_bits(), y.y.to_bits());
                assert_eq!(x.z.to_bits(), y.z.to_bits());
            }
        }
    }

    #[test]
    fn restore_rejects_mismatched_simulation() {
        let mut a = build(false);
        a.run(2, &mut []).unwrap();
        let bytes = a.checkpoint_bytes().unwrap();
        let mut other = crate::integrator::tests::lj_gas(4, 0.7, 1);
        maxwell_boltzmann_velocities(&mut other, 0.8, 4, UnitSystem::Reduced);
        // Smaller cutoff so the 4-cell box is large enough for the neighbour list.
        let ff = LennardJones::single_type(1.0, 1.0, 1.0, Modifier::Shift).unwrap();
        let eval = ForceEvaluator::new(Box::new(ff), 0.2).unwrap();
        let mut small = Simulation::new(
            other,
            eval,
            Box::new(VelocityVerlet),
            None,
            0.002,
            UnitSystem::Reduced,
        )
        .unwrap();
        assert!(small.restore_checkpoint(&bytes).is_err());
    }

    /// Fuzz-style robustness: random bytes, every truncation of a valid file and
    /// bit flips must produce errors, never panics.
    #[test]
    fn decoder_never_panics_on_random_or_truncated_input() {
        let mut sim = build(false);
        sim.run(2, &mut []).unwrap();
        let valid = sim.checkpoint_bytes().unwrap();
        for cut in 0..valid.len() {
            assert!(Checkpoint::decode(&valid[..cut]).is_err());
        }
        let mut rng = Rng::seed_from_u64(0x00C0_FFEE);
        for _ in 0..2000 {
            let len = (rng.next_u64() % 400) as usize;
            let junk: Vec<u8> = (0..len).map(|_| rng.next_u64() as u8).collect();
            let _ = Checkpoint::decode(&junk);
            // Same, but with a valid CRC so the parser itself is exercised.
            let mut framed = junk.clone();
            let crc = crc32(&framed);
            framed.extend_from_slice(&crc.to_le_bytes());
            let _ = Checkpoint::decode(&framed);
        }
        // Valid header with a huge atom count and a correct CRC.
        let mut body = Vec::new();
        body.extend_from_slice(MAGIC);
        body.extend_from_slice(&VERSION.to_le_bytes());
        body.extend_from_slice(&u64::MAX.to_le_bytes());
        body.extend_from_slice(&[0u8; 8 + 8 + 48]);
        let crc = crc32(&body);
        body.extend_from_slice(&crc.to_le_bytes());
        assert!(Checkpoint::decode(&body).is_err());
    }

    #[test]
    fn decoder_rejects_inconsistent_body_lengths() {
        let mut sim = build(false);
        sim.run(1, &mut []).unwrap();
        let valid = sim.checkpoint_bytes().unwrap();
        let body = &valid[..valid.len() - 4];
        // Drop the last state word count by appending trailing garbage.
        let mut extended = body.to_vec();
        extended.push(0);
        assert!(Checkpoint::decode_body(&extended).is_err());
        assert!(Checkpoint::decode_body(&body[..body.len() - 3]).is_err());
    }
}

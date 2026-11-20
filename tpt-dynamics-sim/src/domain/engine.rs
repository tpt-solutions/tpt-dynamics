//! Domain-decomposed Lennard-Jones dynamics with halo exchange, atom migration,
//! and quantile-based load balancing.
//!
//! # Execution model
//!
//! The box is split into `P = nx·ny·nz` sub-domains ([`DomainGrid`]). Each
//! sub-domain is served by one worker thread that owns the atoms whose
//! wrapped positions lie in its region. Workers communicate only through
//! typed `std::sync::mpsc` channels, one channel per ordered pair of workers
//! for each message kind:
//!
//! - **Migration**: after each drift, an atom that now lies in another region
//!   is sent to its new owner.
//! - **Halo (ghost) exchange**: each worker sends every other worker copies of
//!   its atoms that lie within `cutoff` of that worker's region, shifted by the
//!   periodic image that is closest to the region. Ghosts are refreshed every
//!   step, so the neighbour-list skin is not needed.
//!
//! Every send happens before any receive and all channels are unbounded, so the
//! exchange cannot deadlock. A worker that fails drops its senders, which turns
//! the peers' pending receives into errors instead of hangs.
//!
//! Each worker computes forces only on its own atoms and does not use
//! Newton's third law across domains. Every pair within the cutoff is therefore
//! visited from both ends, and each visit carries half the pair energy and
//! virial. This keeps the force kernel free of cross-worker writes.
//!
//! Load balancing runs between segments of `rebalance_every` steps. The
//! coordinator gathers the atoms, recomputes quantile boundaries from the
//! current positions, and scatters the atoms to their new owners. Migration
//! is the only mechanism that moves atoms, so a rebalance is a change of
//! boundaries followed by an ordinary scatter.
//!
//! Rayon is used for the coordinator-side bulk work (the quantile sorts). The
//! workers are plain threads from `std::thread::scope`, because workers that
//! block on channels inside a fixed-size rayon pool can stall if the pool
//! schedules them unevenly.
//!
//! # Scope
//!
//! Supported: orthorhombic boxes, Lennard-Jones with any [`Modifier`],
//! velocity-Verlet NVE. Not supported: thermostats, barostats, constraints,
//! charges and bonded terms. Those belong to the single-domain engine.
//!
//! # Determinism
//!
//! For a fixed grid and schedule, runs are deterministic. Each worker's local
//! ordering is fixed by its sends and receives, and receives are taken in
//! ascending rank order. Forces are not bitwise equal to the single-domain
//! engine, because the summation order differs. Agreement is at the level of
//! floating-point rounding (see the tests).
//!
//! [`Modifier`]: crate::forcefield::Modifier

use std::sync::mpsc::{channel, Receiver, Sender};
use std::time::Instant;

use tpt_dynamics_core::Vec3;

use super::grid::{auto_dims, imbalance, DomainGrid};
use crate::error::SimError;
use crate::forcefield::{ForceField, ForceOutput, LennardJones};
use crate::pbc::SimBox;
use crate::system::System;

/// Options for a decomposed run.
#[derive(Clone, Debug, PartialEq)]
pub struct DecompOptions {
    /// Sub-domain grid. Its product is the number of worker threads.
    pub dims: [usize; 3],
    /// Rebalance every this many steps (0 disables rebalancing).
    pub rebalance_every: usize,
}

impl Default for DecompOptions {
    fn default() -> Self {
        Self { dims: [1, 1, 1], rebalance_every: 0 }
    }
}

impl DecompOptions {
    /// Options for `threads` workers on a box with edges `lengths`, using the
    /// surface-minimising grid from [`auto_dims`].
    pub fn for_threads(threads: usize, lengths: Vec3, rebalance_every: usize) -> Self {
        Self { dims: auto_dims(threads.max(1), lengths), rebalance_every }
    }
}

/// Statistics and final energies of a decomposed run.
#[derive(Clone, Debug, PartialEq)]
pub struct DecompSummary {
    /// Grid used at the start.
    pub dims: [usize; 3],
    /// Steps integrated.
    pub steps: usize,
    /// Potential energy at the final positions.
    pub potential: f64,
    /// Pair virial at the final positions.
    pub virial: f64,
    /// Number of rebalancing operations performed.
    pub rebalances: usize,
    /// Atoms that changed owner, summed over all steps and workers.
    pub migrated: u64,
    /// Ghost (halo) atom copies sent, summed over all steps and workers.
    pub ghosts_sent: u64,
    /// Load imbalance `max/mean` of owned atoms at the end of the run.
    pub imbalance: f64,
    /// Wall-clock time of the run in seconds.
    pub wall_seconds: f64,
    /// Per-phase time, the maximum over workers (critical path), in seconds:
    /// `[migration, ghost exchange, force kernel, kick and drift]`.
    pub phase_seconds: [f64; 4],
}

/// An atom owned by a worker.
#[derive(Clone, Copy, Debug)]
struct Atom {
    id: usize,
    pos: Vec3,
    vel: Vec3,
    force: Vec3,
    mass: f64,
    ty: u16,
}

/// A read-only copy of an atom received from another worker.
#[derive(Clone, Copy, Debug)]
struct Ghost {
    pos: Vec3,
    ty: u16,
}

/// Per-worker channel endpoints. Index `k` is the peer rank.
struct Links {
    mig_tx: Vec<Sender<Vec<Atom>>>,
    mig_rx: Vec<Receiver<Vec<Atom>>>,
    ghost_tx: Vec<Sender<Vec<Ghost>>>,
    ghost_rx: Vec<Receiver<Vec<Ghost>>>,
}

/// Read-only data shared by all workers.
struct Shared<'a> {
    grid: &'a DomainGrid,
    lengths: Vec3,
    cutoff: f64,
    ff: &'a LennardJones,
}

/// Work for one worker over one segment.
struct Job {
    owned: Vec<Atom>,
    steps: usize,
    dt: f64,
    need_eval: bool,
}

/// What a worker returns at the end of a segment.
struct JobResult {
    owned: Vec<Atom>,
    out: ForceOutput,
    migrated: u64,
    ghosts_sent: u64,
    phase: [f64; 4],
}

fn closed() -> SimError {
    SimError::Invalid("domain exchange channel closed: a worker failed".into())
}

fn validate(system: &System, ff: &LennardJones, dt: f64) -> Result<(), SimError> {
    let n = system.len();
    if n == 0 {
        return Err(SimError::Invalid("system has no atoms".into()));
    }
    if system.velocities.len() != n || system.masses.len() != n || system.atom_types.len() != n {
        return Err(SimError::Invalid("per-atom arrays have inconsistent lengths".into()));
    }
    if !(dt.is_finite() && dt > 0.0) {
        return Err(SimError::Invalid(format!("time step must be > 0, got {dt}")));
    }
    if !system.sim_box.is_orthorhombic() {
        return Err(SimError::Invalid("domain decomposition requires an orthorhombic box".into()));
    }
    if system.masses.iter().any(|&m| !(m.is_finite() && m > 0.0)) {
        return Err(SimError::Invalid("all masses must be positive".into()));
    }
    ff.validate_types(&system.atom_types)?;
    let rc = ff.cutoff();
    if 2.0 * rc >= system.sim_box.min_width() {
        return Err(SimError::Invalid(format!(
            "cutoff ({rc} nm) must be less than half the smallest box width ({} nm)",
            system.sim_box.min_width()
        )));
    }
    Ok(())
}

/// Copies the system into atoms with wrapped positions.
fn atoms_from_system(system: &System) -> Vec<Atom> {
    let sim_box: SimBox = system.sim_box;
    (0..system.len())
        .map(|i| Atom {
            id: i,
            pos: sim_box.wrap(system.positions[i]),
            vel: system.velocities[i],
            force: Vec3::ZERO,
            mass: system.masses[i],
            ty: system.atom_types[i],
        })
        .collect()
}

/// Splits atoms into per-domain lists by their current owner.
fn distribute(grid: &DomainGrid, atoms: Vec<Atom>) -> Vec<Vec<Atom>> {
    let mut out: Vec<Vec<Atom>> = (0..grid.n_domains()).map(|_| Vec::new()).collect();
    for a in atoms {
        out[grid.owner(a.pos)].push(a);
    }
    out
}

/// Creates the `P × P` channel mesh. `links[s].mig_tx[d]` and `links[d].mig_rx[s]`
/// are the two ends of the same channel.
fn make_links(p: usize) -> Vec<Links> {
    let mut links: Vec<Links> = (0..p)
        .map(|_| Links {
            mig_tx: Vec::with_capacity(p),
            mig_rx: Vec::with_capacity(p),
            ghost_tx: Vec::with_capacity(p),
            ghost_rx: Vec::with_capacity(p),
        })
        .collect();
    for s in 0..p {
        for d in 0..p {
            let (tx, rx) = channel::<Vec<Atom>>();
            links[s].mig_tx.push(tx);
            links[d].mig_rx.push(rx);
            let (gtx, grx) = channel::<Vec<Ghost>>();
            links[s].ghost_tx.push(gtx);
            links[d].ghost_rx.push(grx);
        }
    }
    links
}

/// Runs one segment on all workers and returns the results in rank order.
fn run_segment(
    grid: &DomainGrid,
    lengths: Vec3,
    ff: &LennardJones,
    owned: Vec<Vec<Atom>>,
    steps: usize,
    dt: f64,
    need_eval: bool,
) -> Result<Vec<JobResult>, SimError> {
    let shared = Shared { grid, lengths, cutoff: ff.cutoff(), ff };
    let links = make_links(grid.n_domains());
    let jobs: Vec<Job> =
        owned.into_iter().map(|owned| Job { owned, steps, dt, need_eval }).collect();
    let shared_ref = &shared;
    std::thread::scope(|scope| {
        let handles: Vec<_> = jobs
            .into_iter()
            .zip(links)
            .enumerate()
            .map(|(rank, (job, lk))| scope.spawn(move || run_worker(rank, shared_ref, job, lk)))
            .collect();
        handles
            .into_iter()
            .map(|h| {
                h.join().unwrap_or_else(|_| Err(SimError::Invalid("domain worker panicked".into())))
            })
            .collect()
    })
}

/// Velocity-Verlet loop for one worker. Every evaluation exchanges migrants and
/// ghosts, so the atoms and their forces are always consistent with the region.
fn run_worker(
    rank: usize,
    shared: &Shared<'_>,
    job: Job,
    links: Links,
) -> Result<JobResult, SimError> {
    let mut w = Worker {
        rank,
        shared,
        links,
        owned: job.owned,
        ghosts: Vec::new(),
        migrated: 0,
        ghosts_sent: 0,
        phase: [0.0; 4],
    };
    let mut out = ForceOutput::default();
    if job.need_eval {
        out = w.evaluate()?;
    }
    let half = 0.5 * job.dt;
    let sim_box = SimBox::new(shared.lengths.x, shared.lengths.y, shared.lengths.z)?;
    for _ in 0..job.steps {
        let t = Instant::now();
        kick(&mut w.owned, half);
        for a in &mut w.owned {
            a.pos = sim_box.wrap(a.pos + a.vel * job.dt);
        }
        w.phase[3] += t.elapsed().as_secs_f64();
        out = w.evaluate()?;
        let t = Instant::now();
        kick(&mut w.owned, half);
        w.phase[3] += t.elapsed().as_secs_f64();
    }
    Ok(JobResult {
        owned: w.owned,
        out,
        migrated: w.migrated,
        ghosts_sent: w.ghosts_sent,
        phase: w.phase,
    })
}

fn kick(atoms: &mut [Atom], factor: f64) {
    for a in atoms {
        a.vel += a.force * (factor / a.mass);
    }
}

struct Worker<'a, 'b> {
    rank: usize,
    shared: &'a Shared<'b>,
    links: Links,
    owned: Vec<Atom>,
    ghosts: Vec<Ghost>,
    migrated: u64,
    ghosts_sent: u64,
    /// Seconds spent in migration, ghost exchange, force kernel, kick and drift.
    phase: [f64; 4],
}

impl Worker<'_, '_> {
    fn evaluate(&mut self) -> Result<ForceOutput, SimError> {
        let t = Instant::now();
        self.migrate()?;
        self.phase[0] += t.elapsed().as_secs_f64();
        let t = Instant::now();
        self.exchange_ghosts()?;
        self.phase[1] += t.elapsed().as_secs_f64();
        let t = Instant::now();
        let out = self.forces();
        self.phase[2] += t.elapsed().as_secs_f64();
        Ok(out)
    }

    /// Sends atoms that left this region to their new owners and receives
    /// atoms that arrived. Every peer gets a message, possibly empty.
    fn migrate(&mut self) -> Result<(), SimError> {
        let grid = self.shared.grid;
        let p = grid.n_domains();
        let mut outgoing: Vec<Vec<Atom>> = (0..p).map(|_| Vec::new()).collect();
        let mut keep = Vec::with_capacity(self.owned.len());
        for a in self.owned.drain(..) {
            let o = grid.owner(a.pos);
            if o == self.rank {
                keep.push(a);
            } else {
                outgoing[o].push(a);
            }
        }
        self.owned = keep;
        for (d, batch) in outgoing.into_iter().enumerate() {
            if d == self.rank {
                continue;
            }
            self.migrated += batch.len() as u64;
            self.links.mig_tx[d].send(batch).map_err(|_| closed())?;
        }
        for s in 0..p {
            if s == self.rank {
                continue;
            }
            let batch = self.links.mig_rx[s].recv().map_err(|_| closed())?;
            self.owned.extend(batch);
        }
        Ok(())
    }

    /// Builds the ghost list: for every target region, the periodic images of
    /// owned atoms that lie within the cutoff of that region.
    fn exchange_ghosts(&mut self) -> Result<(), SimError> {
        let grid = self.shared.grid;
        let lengths = self.shared.lengths;
        let rc = self.shared.cutoff;
        let rc2 = rc * rc;
        let p = grid.n_domains();
        self.ghosts.clear();
        let mut outbound: Vec<Vec<Ghost>> = (0..p).map(|_| Vec::new()).collect();
        let owned = &self.owned;
        let ghosts = &mut self.ghosts;
        let self_rank = self.rank;
        for a in owned {
            for (d, out) in outbound.iter_mut().enumerate() {
                let (lo, hi) = grid.region(d);
                let xs = axis_images(a.pos.x, lo[0], hi[0], lengths.x, rc);
                let ys = axis_images(a.pos.y, lo[1], hi[1], lengths.y, rc);
                let zs = axis_images(a.pos.z, lo[2], hi[2], lengths.z, rc);
                for ix in 0..xs.n {
                    for iy in 0..ys.n {
                        for iz in 0..zs.n {
                            let q = Vec3::new(xs.val[ix], ys.val[iy], zs.val[iz]);
                            if box_dist2(q, lo, hi) >= rc2 {
                                continue;
                            }
                            let g = Ghost { pos: q, ty: a.ty };
                            if d == self_rank {
                                // The unshifted atom itself is not its own ghost.
                                if xs.zero[ix] && ys.zero[iy] && zs.zero[iz] {
                                    continue;
                                }
                                ghosts.push(g);
                            } else {
                                out.push(g);
                            }
                        }
                    }
                }
            }
        }
        for (d, batch) in outbound.into_iter().enumerate() {
            if d == self.rank {
                continue;
            }
            self.ghosts_sent += batch.len() as u64;
            self.links.ghost_tx[d].send(batch).map_err(|_| closed())?;
        }
        for s in 0..p {
            if s == self.rank {
                continue;
            }
            let batch = self.links.ghost_rx[s].recv().map_err(|_| closed())?;
            self.ghosts.extend(batch);
        }
        Ok(())
    }

    /// Lennard-Jones forces on owned atoms from owned and ghost atoms, using a
    /// local cell grid. Each pair is visited from both ends with half the energy.
    fn forces(&mut self) -> ForceOutput {
        let ff = self.shared.ff;
        let rc = self.shared.cutoff;
        let rc2 = rc * rc;
        let n_own = self.owned.len();
        let mut out = ForceOutput::default();
        if n_own == 0 {
            return out;
        }
        let mut pos: Vec<Vec3> = Vec::with_capacity(n_own + self.ghosts.len());
        let mut ty: Vec<u16> = Vec::with_capacity(n_own + self.ghosts.len());
        for a in &self.owned {
            pos.push(a.pos);
            ty.push(a.ty);
        }
        for g in &self.ghosts {
            pos.push(g.pos);
            ty.push(g.ty);
        }

        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        for p in &pos {
            for (a, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                lo[a] = lo[a].min(v);
                hi[a] = hi[a].max(v);
            }
        }
        // Cells at least `rc` wide: the 27 neighbouring cells cover the cutoff sphere.
        let mut nc = [1_usize; 3];
        let mut cs = [1.0_f64; 3];
        for a in 0..3 {
            let ext = hi[a] - lo[a];
            nc[a] = ((ext / rc).floor() as usize).max(1);
            cs[a] = if ext > 0.0 { ext / nc[a] as f64 } else { 1.0 };
        }
        let cell_id = |c: [usize; 3]| c[0] + nc[0] * (c[1] + nc[1] * c[2]);
        let cell_xyz: Vec<[usize; 3]> = pos
            .iter()
            .map(|p| {
                let mut c = [0_usize; 3];
                for (a, v) in [p.x, p.y, p.z].into_iter().enumerate() {
                    let k = ((v - lo[a]) / cs[a]).floor().max(0.0) as usize;
                    c[a] = k.min(nc[a] - 1);
                }
                c
            })
            .collect();
        let ncell = nc[0] * nc[1] * nc[2];
        let mut start = vec![0_usize; ncell + 1];
        for c in &cell_xyz {
            start[cell_id(*c) + 1] += 1;
        }
        for k in 0..ncell {
            start[k + 1] += start[k];
        }
        let mut fill = start.clone();
        let mut members = vec![0_usize; pos.len()];
        for (j, c) in cell_xyz.iter().enumerate() {
            let id = cell_id(*c);
            members[fill[id]] = j;
            fill[id] += 1;
        }

        let mut forces = vec![Vec3::ZERO; n_own];
        for i in 0..n_own {
            let c = cell_xyz[i];
            let x_range = c[0].saturating_sub(1)..=(c[0] + 1).min(nc[0] - 1);
            let y_range = c[1].saturating_sub(1)..=(c[1] + 1).min(nc[1] - 1);
            let z_range = c[2].saturating_sub(1)..=(c[2] + 1).min(nc[2] - 1);
            let mut fi = Vec3::ZERO;
            for cz in z_range {
                for cy in y_range.clone() {
                    for cx in x_range.clone() {
                        let id = cell_id([cx, cy, cz]);
                        for &j in &members[start[id]..start[id + 1]] {
                            if j == i {
                                continue;
                            }
                            let dr = pos[i] - pos[j];
                            let r2 = dr.norm_sq();
                            if r2 >= rc2 || r2 == 0.0 {
                                continue;
                            }
                            let r = r2.sqrt();
                            let (v, dv) = ff.pair(ty[i], ty[j], r);
                            let f = dr * (-dv / r);
                            fi += f;
                            out.potential += 0.5 * v;
                            out.virial += 0.5 * dr.dot(f);
                        }
                    }
                }
            }
            forces[i] = fi;
        }
        for (a, f) in self.owned.iter_mut().zip(forces) {
            a.force = f;
        }
        out
    }
}

/// Periodic images of coordinate `x` (in `[0, len)`) that lie within `rc` of
/// the interval `[lo, hi]`. Shifts of `-1, 0, +1` box lengths suffice when
/// `2·rc < len`.
struct AxisImages {
    val: [f64; 3],
    zero: [bool; 3],
    n: usize,
}

fn axis_images(x: f64, lo: f64, hi: f64, len: f64, rc: f64) -> AxisImages {
    let mut out = AxisImages { val: [0.0; 3], zero: [false; 3], n: 0 };
    for s in [-1.0_f64, 0.0, 1.0] {
        let c = x + s * len;
        let d = if c < lo {
            lo - c
        } else if c > hi {
            c - hi
        } else {
            0.0
        };
        if d < rc {
            out.val[out.n] = c;
            out.zero[out.n] = s == 0.0;
            out.n += 1;
        }
    }
    out
}

/// Squared distance from `q` to the box `[lo, hi]` (zero if inside).
fn box_dist2(q: Vec3, lo: [f64; 3], hi: [f64; 3]) -> f64 {
    let mut d2 = 0.0;
    for (a, c) in [q.x, q.y, q.z].into_iter().enumerate() {
        let d = if c < lo[a] {
            lo[a] - c
        } else if c > hi[a] {
            c - hi[a]
        } else {
            0.0
        };
        d2 += d * d;
    }
    d2
}

/// Computes potential, virial and forces for the current positions of `system`
/// with a single decomposed evaluation over `dims` domains. Does not change
/// `system`. Forces are returned in atom order.
///
/// # Errors
/// Returns an error for systems the decomposition does not support (see the
/// module docs) or if the box is too small for the cutoff.
pub fn decomposed_forces(
    system: &System,
    ff: &LennardJones,
    dims: [usize; 3],
) -> Result<(ForceOutput, Vec<Vec3>), SimError> {
    validate(system, ff, 1.0)?;
    let lengths = system.sim_box.lengths;
    let grid = DomainGrid::uniform(dims, lengths);
    let owned = distribute(&grid, atoms_from_system(system));
    let results = run_segment(&grid, lengths, ff, owned, 0, 0.0, true)?;
    let mut forces = vec![Vec3::ZERO; system.len()];
    let mut total = ForceOutput::default();
    for r in results {
        total.potential += r.out.potential;
        total.virial += r.out.virial;
        for a in r.owned {
            forces[a.id] = a.force;
        }
    }
    Ok((total, forces))
}

/// Integrates `nsteps` of velocity-Verlet NVE dynamics with the decomposed
/// engine. Positions are wrapped into the box at the start. Velocities and
/// forces are written back to `system`.
///
/// # Errors
/// Returns an error for unsupported systems, an invalid time step, or a
/// failed exchange.
pub fn run_decomposed_nve(
    system: &mut System,
    ff: &LennardJones,
    dt: f64,
    nsteps: usize,
    opts: &DecompOptions,
) -> Result<DecompSummary, SimError> {
    validate(system, ff, dt)?;
    let t0 = Instant::now();
    let lengths = system.sim_box.lengths;
    let mut grid = DomainGrid::uniform(opts.dims, lengths);
    let mut atoms = atoms_from_system(system);
    let mut owned = distribute(&grid, atoms.clone());
    let mut remaining = nsteps;
    let mut need_eval = true;
    let mut rebalances = 0;
    let mut migrated = 0_u64;
    let mut ghosts_sent = 0_u64;
    let mut phase = [0.0_f64; 4];
    let mut last;
    loop {
        let seg =
            if opts.rebalance_every == 0 { remaining } else { remaining.min(opts.rebalance_every) };
        let results = run_segment(&grid, lengths, ff, owned, seg, dt, need_eval)?;
        need_eval = false;
        remaining -= seg;
        let mut seg_out = ForceOutput::default();
        let mut by_id: Vec<Option<Atom>> = vec![None; atoms.len()];
        for r in results {
            seg_out.potential += r.out.potential;
            seg_out.virial += r.out.virial;
            migrated += r.migrated;
            ghosts_sent += r.ghosts_sent;
            for (acc, t) in phase.iter_mut().zip(r.phase) {
                *acc = acc.max(t);
            }
            for a in r.owned {
                by_id[a.id] = Some(a);
            }
        }
        last = seg_out;
        atoms = by_id
            .into_iter()
            .map(|a| a.ok_or_else(|| SimError::Invalid("atom lost during exchange".into())))
            .collect::<Result<_, _>>()?;
        if remaining == 0 {
            break;
        }
        if opts.rebalance_every > 0 {
            let positions: Vec<Vec3> = atoms.iter().map(|a| a.pos).collect();
            grid = grid.rebalanced(&positions, lengths);
            rebalances += 1;
        }
        owned = distribute(&grid, atoms.clone());
    }
    let mut counts = vec![0_usize; grid.n_domains()];
    for a in &atoms {
        counts[grid.owner(a.pos)] += 1;
    }
    // Write back. Atoms are in id order, so the arrays stay aligned with `system`.
    for a in &atoms {
        system.positions[a.id] = a.pos;
        system.velocities[a.id] = a.vel;
        system.forces[a.id] = a.force;
    }
    Ok(DecompSummary {
        dims: opts.dims,
        steps: nsteps,
        potential: last.potential,
        virial: last.virial,
        rebalances,
        migrated,
        ghosts_sent,
        imbalance: imbalance(&counts),
        wall_seconds: t0.elapsed().as_secs_f64(),
        phase_seconds: phase,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::Modifier;
    use crate::integrator::{ForceEvaluator, VelocityVerlet};
    use crate::rng::Rng;
    use crate::simulation::Simulation;
    use crate::system::{NewAtom, UnitSystem};

    /// FCC crystal of `cells^3` conventional cells at reduced density 0.8442,
    /// with a seeded random displacement and velocity so that atoms migrate.
    fn fcc_system(cells: usize, seed: u64) -> System {
        let a = (4.0_f64 / 0.8442).cbrt();
        let l = cells as f64 * a;
        let mut sys = System::new(SimBox::new(l, l, l).unwrap());
        let basis = [
            Vec3::new(0.0, 0.0, 0.0),
            Vec3::new(0.5, 0.5, 0.0),
            Vec3::new(0.5, 0.0, 0.5),
            Vec3::new(0.0, 0.5, 0.5),
        ];
        let mut rng = Rng::seed_from_u64(seed);
        for ix in 0..cells {
            for iy in 0..cells {
                for iz in 0..cells {
                    for b in basis {
                        let p = Vec3::new(ix as f64 + b.x, iy as f64 + b.y, iz as f64 + b.z) * a;
                        let jitter = Vec3::new(
                            0.05 * (rng.uniform() - 0.5),
                            0.05 * (rng.uniform() - 0.5),
                            0.05 * (rng.uniform() - 0.5),
                        );
                        let v = Vec3::new(
                            rng.uniform() - 0.5,
                            rng.uniform() - 0.5,
                            rng.uniform() - 0.5,
                        );
                        sys.push_atom(NewAtom {
                            name: "AR".into(),
                            residue_name: "LJ".into(),
                            residue_id: 1,
                            element: "Ar".into(),
                            position: p + jitter,
                            velocity: v,
                            mass: 1.0,
                            charge: 0.0,
                            atom_type: 0,
                        });
                    }
                }
            }
        }
        sys
    }

    fn lj() -> LennardJones {
        LennardJones::single_type(1.0, 1.0, 2.5, Modifier::Shift).unwrap()
    }

    fn min_image_dist(sb: &SimBox, a: Vec3, b: Vec3) -> f64 {
        sb.min_image(a - b).norm()
    }

    #[test]
    fn forces_and_energy_match_single_domain() {
        let sys = fcc_system(4, 7);
        let ff = lj();
        let mut reference = sys.clone();
        let mut eval = ForceEvaluator::new(Box::new(ff.clone()), 0.3).unwrap();
        let single = eval.evaluate(&mut reference).unwrap();
        for dims in [[1, 1, 1], [2, 1, 1], [2, 2, 2], [3, 2, 1], [4, 1, 1], [1, 4, 2]] {
            let (out, forces) = decomposed_forces(&sys, &ff, dims).unwrap();
            let de = (out.potential - single.potential).abs();
            assert!(
                de <= 1e-10 * single.potential.abs(),
                "{dims:?}: potential {} vs {}",
                out.potential,
                single.potential
            );
            assert!((out.virial - single.virial).abs() <= 1e-10 * single.virial.abs().max(1.0));
            for (i, (f, g)) in forces.iter().zip(&reference.forces).enumerate() {
                let diff = (*f - *g).norm();
                assert!(diff <= 1e-9 * (1.0 + g.norm()), "{dims:?} atom {i}: {diff}");
            }
        }
    }

    #[test]
    fn nve_trajectory_matches_single_domain_with_rebalancing() {
        let sys = fcc_system(4, 11);
        let ff = lj();
        let dt = 0.002;
        let steps = 20;
        let eval = ForceEvaluator::new(Box::new(ff.clone()), 0.3).unwrap();
        let mut sim = Simulation::new(
            sys.clone(),
            eval,
            Box::new(VelocityVerlet),
            None,
            dt,
            UnitSystem::Reduced,
        )
        .unwrap();
        for _ in 0..steps {
            sim.advance().unwrap();
        }
        let mut dec = sys.clone();
        let opts = DecompOptions { dims: [2, 2, 1], rebalance_every: 5 };
        let summary = run_decomposed_nve(&mut dec, &ff, dt, steps, &opts).unwrap();
        assert_eq!(summary.rebalances, 3, "rebalance after steps 5, 10, 15");
        assert!(summary.migrated > 0, "atoms should have crossed domain boundaries");
        assert!(summary.ghosts_sent > 0);
        for i in 0..sys.len() {
            let d = min_image_dist(&sys.sim_box, dec.positions[i], sim.system.positions[i]);
            assert!(d < 1e-9, "atom {i}: position differs by {d}");
            let dv = (dec.velocities[i] - sim.system.velocities[i]).norm();
            assert!(dv < 1e-9, "atom {i}: velocity differs by {dv}");
        }
        let e_ref = sim.energies();
        assert!((summary.potential - e_ref.potential).abs() <= 1e-9 * e_ref.potential.abs());
    }

    #[test]
    fn repeated_decomposed_runs_are_bitwise_identical() {
        let sys = fcc_system(4, 3);
        let ff = lj();
        let opts = DecompOptions { dims: [2, 2, 2], rebalance_every: 4 };
        let mut a = sys.clone();
        let mut b = sys.clone();
        let sa = run_decomposed_nve(&mut a, &ff, 0.002, 12, &opts).unwrap();
        let sb = run_decomposed_nve(&mut b, &ff, 0.002, 12, &opts).unwrap();
        assert_eq!(a.positions, b.positions);
        assert_eq!(a.velocities, b.velocities);
        assert_eq!(sa.potential.to_bits(), sb.potential.to_bits());
        assert_eq!(sa.migrated, sb.migrated);
    }

    #[test]
    fn single_worker_has_no_exchange_traffic() {
        let sys = fcc_system(4, 5);
        let (_, _) = decomposed_forces(&sys, &lj(), [1, 1, 1]).unwrap();
        let mut s = sys.clone();
        let r = run_decomposed_nve(&mut s, &lj(), 0.002, 3, &DecompOptions::default()).unwrap();
        assert_eq!(r.migrated, 0);
        assert_eq!(r.rebalances, 0);
    }

    #[test]
    fn rejects_box_too_small_for_cutoff() {
        let sys = fcc_system(2, 1);
        let ff = LennardJones::single_type(1.0, 1.0, 2.5, Modifier::None).unwrap();
        assert!(decomposed_forces(&sys, &ff, [1, 1, 1]).is_err());
    }
}

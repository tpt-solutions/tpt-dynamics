//! Simple kernel timing harness (no external benchmark framework).
//!
//! Run with `cargo run --release -p tpt-dynamics-compute --example kernel_bench`.
//! Each kernel is run a fixed number of times and the minimum and median wall-clock
//! times are reported. Numbers depend heavily on the CPU and build flags; see
//! `docs/benchmarks/compute-kernels.md` for the recorded results and the machine.

#![allow(clippy::many_single_char_names, clippy::similar_names)]

use std::time::Instant;

use tpt_dynamics_compute::cluster::{build_cluster_pairs, cluster_pair_kernel, ClusterSet};
use tpt_dynamics_compute::pme::{Pme, PmeParams};
use tpt_dynamics_compute::scalar::pair_list_scalar;
use tpt_dynamics_compute::simd_kernel::pair_list_simd;
use tpt_dynamics_compute::{AtomParams, ForceBuffers, NonbondedConfig, Positions};

/// Simple xorshift generator in [0, 1).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

/// Times `f` over `reps` runs; returns (min, median) in seconds.
fn time<F: FnMut()>(reps: usize, mut f: F) -> (f64, f64) {
    f(); // warm-up
    let mut t: Vec<f64> = (0..reps)
        .map(|_| {
            let s = Instant::now();
            f();
            s.elapsed().as_secs_f64()
        })
        .collect();
    t.sort_by(|a, b| a.partial_cmp(b).unwrap());
    (t[0], t[t.len() / 2])
}

fn random_positions(n: usize, edge: f64, seed: u64) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut rng = Rng(seed);
    let mut v = (Vec::with_capacity(n), Vec::with_capacity(n), Vec::with_capacity(n));
    for _ in 0..n {
        v.0.push(rng.next() * edge);
        v.1.push(rng.next() * edge);
        v.2.push(rng.next() * edge);
    }
    v
}

#[allow(clippy::too_many_lines)]
fn main() {
    let reps = 9;
    let threads = rayon::current_num_threads();
    println!("threads (rayon): {threads}");

    // Pair-kernel system: LJ-fluid density (0.8 reduced), cutoff 2.5.
    let n = 4096;
    let rc = 2.5;
    let edge = (n as f64 / 0.8).cbrt();
    let (x, y, z) = random_positions(n, edge, 2024);
    let sigma = vec![1.0; n];
    let eps = vec![1.0; n];
    let q: Vec<f64> = (0..n).map(|i| if i % 2 == 0 { 1.0 } else { -1.0 }).collect();
    let pos = Positions::new(&x, &y, &z);
    let params = AtomParams { sigma: &sigma, epsilon: &eps, charge: &q };
    let cfg = NonbondedConfig { cutoff: rc, alpha: 1.4, box_len: Some([edge; 3]) };

    // Verlet-style pair list built once (brute force, outside the timed region).
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            let (dx, dy, dz) = (x[i] - x[j], y[i] - y[j], z[i] - z[j]);
            let wrap = |d: f64| d - edge * (d / edge).round();
            let r2 = wrap(dx).powi(2) + wrap(dy).powi(2) + wrap(dz).powi(2);
            if r2 < rc * rc {
                pairs.push((i as u32, j as u32));
            }
        }
    }
    println!("system: N={n}, box={edge:.3}, cutoff={rc}, pairs={}", pairs.len());

    let mut fx = vec![0.0; n];
    let mut fy = vec![0.0; n];
    let mut fz = vec![0.0; n];
    let (smin, smed) = time(reps, || {
        fx.fill(0.0);
        fy.fill(0.0);
        fz.fill(0.0);
        let _ = pair_list_scalar(
            pos,
            params,
            &pairs,
            &cfg,
            &mut ForceBuffers::new(&mut fx, &mut fy, &mut fz),
        );
    });
    report("pair list, scalar", smin, smed, pairs.len());

    let (vmin, vmed) = time(reps, || {
        fx.fill(0.0);
        fy.fill(0.0);
        fz.fill(0.0);
        let _ = pair_list_simd(
            pos,
            params,
            &pairs,
            &cfg,
            &mut ForceBuffers::new(&mut fx, &mut fy, &mut fz),
        );
    });
    report("pair list, tpt-simd", vmin, vmed, pairs.len());

    let set = ClusterSet::from_atoms(pos, params);
    let cpairs = build_cluster_pairs(&set, rc, cfg.box_len);
    let mut cf = set.zero_forces();
    let (cmin, cmed) = time(reps, || {
        for v in cf.fx.iter_mut().chain(cf.fy.iter_mut()).chain(cf.fz.iter_mut()) {
            *v = [0.0; 4];
        }
        let _ = cluster_pair_kernel(&set, &cpairs, &cfg, &mut cf);
    });
    report(&format!("cluster pairs ({} cluster pairs)", cpairs.len()), cmin, cmed, pairs.len());

    // Same cluster kernel after a spatial sort (cells of size rc), as a simulation would do.
    let mut order: Vec<usize> = (0..n).collect();
    let cell = |v: f64| (v.rem_euclid(edge) / rc).floor() as i64;
    order.sort_by_key(|&a| (cell(x[a]), cell(y[a]), cell(z[a])));
    let gather = |v: &[f64]| order.iter().map(|&a| v[a]).collect::<Vec<f64>>();
    let (sx, sy, sz) = (gather(&x), gather(&y), gather(&z));
    let (ss, se, sq) = (gather(&sigma), gather(&eps), gather(&q));
    let sset = ClusterSet::from_atoms(
        Positions::new(&sx, &sy, &sz),
        AtomParams { sigma: &ss, epsilon: &se, charge: &sq },
    );
    let spairs = build_cluster_pairs(&sset, rc, cfg.box_len);
    let mut scf = sset.zero_forces();
    let (tmin, tmed) = time(reps, || {
        for v in scf.fx.iter_mut().chain(scf.fy.iter_mut()).chain(scf.fz.iter_mut()) {
            *v = [0.0; 4];
        }
        let _ = cluster_pair_kernel(&sset, &spairs, &cfg, &mut scf);
    });
    report(
        &format!("cluster pairs, spatially sorted ({} pairs)", spairs.len()),
        tmin,
        tmed,
        pairs.len(),
    );

    // PME reciprocal solve on a smaller box so the grid stays modest in memory.
    let np = 512;
    let pbox = 8.0;
    let (px, py, pz) = random_positions(np, pbox, 99);
    let pq: Vec<f64> = (0..np).map(|i| if i % 2 == 0 { 1.0 } else { -1.0 }).collect();
    let mesh_pos = Positions::new(&px, &py, &pz);
    let pme = Pme::new([pbox; 3], PmeParams { alpha: 2.0, grid: [64; 3], order: 4 });
    let mut gx = vec![0.0; np];
    let mut gy = vec![0.0; np];
    let mut gz = vec![0.0; np];
    let (pmin, pmed) = time(reps, || {
        gx.fill(0.0);
        gy.fill(0.0);
        gz.fill(0.0);
        let _ = pme.reciprocal(mesh_pos, &pq, &mut ForceBuffers::new(&mut gx, &mut gy, &mut gz));
    });
    report(&format!("PME reciprocal (N={np}, grid 64^3, order 4)"), pmin, pmed, np);
}

fn report(name: &str, min: f64, median: f64, items: usize) {
    println!(
        "{name:<48} min {:>9.3} ms  median {:>9.3} ms  {:>8.2} ns/item (median)",
        min * 1e3,
        median * 1e3,
        median * 1e9 / items as f64
    );
}

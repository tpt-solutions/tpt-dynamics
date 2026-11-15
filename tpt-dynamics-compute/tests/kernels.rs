//! Correctness of the pairwise kernels: scalar reference vs SIMD vs cluster layout,
//! and force consistency with the energy formulas.

#![allow(clippy::many_single_char_names, clippy::similar_names, clippy::type_complexity)]

use tpt_dynamics_compute::cluster::{build_cluster_pairs, cluster_pair_kernel, ClusterSet};
use tpt_dynamics_compute::pairwise::{coulomb_ewald_pair, excluded_ewald_pair, lj_pair};
use tpt_dynamics_compute::scalar::{excluded_pairs_correction, pair_list_scalar};
use tpt_dynamics_compute::simd_kernel::pair_list_simd;
use tpt_dynamics_compute::{AtomParams, Energies, ForceBuffers, NonbondedConfig, Positions};

/// Deterministic xorshift generator in [0, 1).
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
}

struct System {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    sigma: Vec<f64>,
    eps: Vec<f64>,
    q: Vec<f64>,
}

fn random_system(n: usize, edge: f64, seed: u64) -> System {
    let mut rng = Rng(seed);
    let mut s = System { x: vec![], y: vec![], z: vec![], sigma: vec![], eps: vec![], q: vec![] };
    for i in 0..n {
        s.x.push(rng.next() * edge);
        s.y.push(rng.next() * edge);
        s.z.push(rng.next() * edge);
        s.sigma.push(0.25 + 0.2 * rng.next());
        s.eps.push(0.5 + rng.next());
        // Alternating charges, slightly perturbed, so the system is not trivially neutral per pair.
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        s.q.push(sign * (0.5 + rng.next()));
    }
    s
}

fn all_pairs(n: usize) -> Vec<(u32, u32)> {
    let mut v = Vec::new();
    for i in 0..n {
        for j in i + 1..n {
            v.push((i as u32, j as u32));
        }
    }
    v
}

fn assert_close(a: f64, b: f64, tol: f64, what: &str) {
    let scale = a.abs().max(b.abs()).max(1.0);
    assert!((a - b).abs() <= tol * scale, "{what}: {a} vs {b}");
}

#[test]
fn lj_minimum_is_minus_epsilon() {
    let sigma = 1.3;
    let rmin = 2f64.powf(1.0 / 6.0) * sigma;
    let (e, f) = lj_pair(rmin * rmin, sigma, 0.7);
    assert!((e + 0.7).abs() < 1e-12, "energy at minimum: {e}");
    assert!(f.abs() < 1e-10, "force at minimum: {f}");
}

#[test]
fn pair_formulas_match_energy_derivative() {
    // Radial force F = -dE/dr, and f_over_r * r = F. Central differences.
    let h = 1e-6;
    let alpha = 2.1;
    for &r in &[0.6_f64, 1.0, 1.7] {
        let r2 = r * r;
        let energy = |rr: f64| {
            let (e1, _) = lj_pair(rr * rr, 0.9, 0.5);
            let (e2, _) = coulomb_ewald_pair(rr * rr, -0.8, alpha);
            let (e3, _) = excluded_ewald_pair(rr * rr, 0.6, alpha);
            (e1, e2, e3)
        };
        let (p, m) = (energy(r + h), energy(r - h));
        let (f_lj, f_c, f_ex) =
            (-(p.0 - m.0) / (2.0 * h), -(p.1 - m.1) / (2.0 * h), -(p.2 - m.2) / (2.0 * h));
        let (_, flj) = lj_pair(r2, 0.9, 0.5);
        let (_, fc) = coulomb_ewald_pair(r2, -0.8, alpha);
        let (_, fe) = excluded_ewald_pair(r2, 0.6, alpha);
        assert!((flj * r - f_lj).abs() < 1e-6, "LJ force at r={r}");
        assert!((fc * r - f_c).abs() < 1e-6, "Coulomb force at r={r}");
        assert!((fe * r - f_ex).abs() < 1e-6, "excluded force at r={r}");
    }
}

#[test]
fn simd_pair_list_matches_scalar() {
    let n = 203; // not a multiple of 4 or of the cluster size
    let edge = 5.0;
    let s = random_system(n, edge, 42);
    let pos = Positions::new(&s.x, &s.y, &s.z);
    let params = AtomParams { sigma: &s.sigma, epsilon: &s.eps, charge: &s.q };
    let pairs = all_pairs(n);
    let cfg = NonbondedConfig { cutoff: 1.6, alpha: 2.4, box_len: Some([edge; 3]) };

    let mut fs = vec![[0.0; 3]; n];
    let mut fv = vec![[0.0; 3]; n];
    let (mut ax, mut ay, mut az) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let (mut bx, mut by, mut bz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    let es = pair_list_scalar(
        pos,
        params,
        &pairs,
        &cfg,
        &mut ForceBuffers::new(&mut ax, &mut ay, &mut az),
    );
    let ev = pair_list_simd(
        pos,
        params,
        &pairs,
        &cfg,
        &mut ForceBuffers::new(&mut bx, &mut by, &mut bz),
    );
    for i in 0..n {
        fs[i] = [ax[i], ay[i], az[i]];
        fv[i] = [bx[i], by[i], bz[i]];
    }
    assert!(es.lj != 0.0 && es.coulomb != 0.0, "test system should have interacting pairs");
    assert_close(es.lj, ev.lj, 1e-12, "LJ energy");
    assert_close(es.coulomb, ev.coulomb, 1e-12, "Coulomb energy");
    for i in 0..n {
        for k in 0..3 {
            assert_close(fs[i][k], fv[i][k], 1e-11, "force");
        }
    }
}

#[test]
fn cluster_kernel_matches_scalar_with_padding() {
    for &n in &[4usize, 13, 37] {
        let edge = 4.0;
        let s = random_system(n, edge, 7 + n as u64);
        let pos = Positions::new(&s.x, &s.y, &s.z);
        let params = AtomParams { sigma: &s.sigma, epsilon: &s.eps, charge: &s.q };
        let cfg = NonbondedConfig { cutoff: 1.4, alpha: 2.0, box_len: Some([edge; 3]) };

        // Scalar reference over all atom pairs.
        let (mut ax, mut ay, mut az) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        let es = pair_list_scalar(
            pos,
            params,
            &all_pairs(n),
            &cfg,
            &mut ForceBuffers::new(&mut ax, &mut ay, &mut az),
        );

        // Cluster path.
        let set = ClusterSet::from_atoms(pos, params);
        let cpairs = build_cluster_pairs(&set, cfg.cutoff, cfg.box_len);
        let mut cf = set.zero_forces();
        let ec: Energies = cluster_pair_kernel(&set, &cpairs, &cfg, &mut cf);
        let (mut bx, mut by, mut bz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
        set.scatter_forces(&cf, &mut ForceBuffers::new(&mut bx, &mut by, &mut bz));

        assert_close(es.lj, ec.lj, 1e-12, &format!("LJ energy n={n}"));
        assert_close(es.coulomb, ec.coulomb, 1e-12, &format!("Coulomb energy n={n}"));
        for i in 0..n {
            assert_close(ax[i], bx[i], 1e-11, "fx");
            assert_close(ay[i], by[i], 1e-11, "fy");
            assert_close(az[i], bz[i], 1e-11, "fz");
        }
    }
}

#[test]
fn forces_obey_newton_third_law() {
    let n = 60;
    let edge = 4.0;
    let s = random_system(n, edge, 99);
    let pos = Positions::new(&s.x, &s.y, &s.z);
    let params = AtomParams { sigma: &s.sigma, epsilon: &s.eps, charge: &s.q };
    let cfg = NonbondedConfig { cutoff: 1.8, alpha: 2.0, box_len: Some([edge; 3]) };
    let (mut fx, mut fy, mut fz) = (vec![0.0; n], vec![0.0; n], vec![0.0; n]);
    pair_list_simd(
        pos,
        params,
        &all_pairs(n),
        &cfg,
        &mut ForceBuffers::new(&mut fx, &mut fy, &mut fz),
    );
    let total: f64 = fx.iter().chain(&fy).chain(&fz).map(|v| v.abs()).sum();
    assert!(total > 0.1);
    for v in [&fx, &fy, &fz] {
        assert!(v.iter().sum::<f64>().abs() < 1e-10, "net force not zero");
    }
}

#[test]
fn excluded_pair_correction_matches_numeric_gradient() {
    // Two atoms, open geometry (no box): check forces against energy differences.
    let x = [0.0, 0.7];
    let y = [0.0, 0.2];
    let z = [0.0, -0.3];
    let q = [1.0, -0.6];
    let pos = Positions::new(&x, &y, &z);
    let alpha = 1.9;
    let (mut fx, mut fy, mut fz) = (vec![0.0; 2], vec![0.0; 2], vec![0.0; 2]);
    let e0 = excluded_pairs_correction(
        pos,
        &q,
        &[(0, 1)],
        alpha,
        None,
        &mut ForceBuffers::new(&mut fx, &mut fy, &mut fz),
    );
    let h = 1e-6;
    let mut xp = x;
    xp[0] += h;
    let mut junk = (vec![0.0; 2], vec![0.0; 2], vec![0.0; 2]);
    let ep = excluded_pairs_correction(
        Positions::new(&xp, &y, &z),
        &q,
        &[(0, 1)],
        alpha,
        None,
        &mut ForceBuffers::new(&mut junk.0, &mut junk.1, &mut junk.2),
    );
    let numeric = -(ep - e0) / h;
    assert!((fx[0] - numeric).abs() < 1e-5, "excluded force {} vs numeric {numeric}", fx[0]);
    assert!((fx[0] + fx[1]).abs() < 1e-12);
}

//! PME accuracy against brute-force Ewald summation on a small neutral system.
//!
//! System: 12 unit-magnitude charges (alternating sign) in a cubic box of edge 2.0,
//! Ewald alpha 3.8 and real-space cutoff 1.0 (= L/2), so erfc(alpha * rc) ~ 7e-8.
//! Reference: [`direct_reciprocal`] and [`direct_real_images`] with converged sums.
//! Tolerances are stated per test; they bound the PME mesh error for grid 32, order 6.

#![allow(clippy::many_single_char_names, clippy::similar_names, clippy::type_complexity)]

use tpt_dynamics_compute::cluster::{build_cluster_pairs, cluster_pair_kernel, ClusterSet};
use tpt_dynamics_compute::pme::{alpha_for_tolerance, Pme, PmeParams};
use tpt_dynamics_compute::reference::{direct_real_images, direct_reciprocal, direct_self_energy};
use tpt_dynamics_compute::scalar::pair_list_scalar;
use tpt_dynamics_compute::simd_kernel::pair_list_simd;
use tpt_dynamics_compute::{AtomParams, ForceBuffers, NonbondedConfig, Positions};

const EDGE: f64 = 2.0;
const N: usize = 12;

struct Neutral {
    x: Vec<f64>,
    y: Vec<f64>,
    z: Vec<f64>,
    q: Vec<f64>,
}

fn neutral_system(seed: u64) -> Neutral {
    let mut s = seed;
    let mut next = move || {
        s ^= s << 13;
        s ^= s >> 7;
        s ^= s << 17;
        (s >> 11) as f64 / (1u64 << 53) as f64
    };
    let mut sys = Neutral { x: vec![], y: vec![], z: vec![], q: vec![] };
    for i in 0..N {
        sys.x.push(next() * EDGE);
        sys.y.push(next() * EDGE);
        sys.z.push(next() * EDGE);
        sys.q.push(if i % 2 == 0 { 1.0 } else { -1.0 });
    }
    sys
}

fn zeros() -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    (vec![0.0; N], vec![0.0; N], vec![0.0; N])
}

fn rel_err(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

/// RMS of force differences relative to RMS of reference forces.
fn rms_rel(a: &(Vec<f64>, Vec<f64>, Vec<f64>), b: &(Vec<f64>, Vec<f64>, Vec<f64>)) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for k in 0..N {
        for (va, vb) in [(&a.0, &b.0), (&a.1, &b.1), (&a.2, &b.2)] {
            num += (va[k] - vb[k]).powi(2);
            den += vb[k].powi(2);
        }
    }
    (num / den).sqrt()
}

fn direct_total(sys: &Neutral, alpha: f64) -> (f64, (Vec<f64>, Vec<f64>, Vec<f64>)) {
    let pos = Positions::new(&sys.x, &sys.y, &sys.z);
    let box_len = [EDGE; 3];
    let mut f = zeros();
    let e_real = direct_real_images(
        pos,
        &sys.q,
        box_len,
        alpha,
        2,
        &mut ForceBuffers::new(&mut f.0, &mut f.1, &mut f.2),
    );
    let e_rec = direct_reciprocal(
        pos,
        &sys.q,
        box_len,
        alpha,
        10,
        &mut ForceBuffers::new(&mut f.0, &mut f.1, &mut f.2),
    );
    (e_real + e_rec + direct_self_energy(&sys.q, alpha), f)
}

#[test]
fn direct_ewald_is_independent_of_alpha() {
    // The total Ewald energy of a neutral periodic system does not depend on alpha.
    // This validates the brute-force reference before it is used to judge PME.
    let sys = neutral_system(1234);
    let (e1, _) = direct_total(&sys, 3.8);
    let (e2, _) = direct_total(&sys, 2.9);
    assert!((e1 - e2).abs() < 1e-6, "alpha dependence: {e1} vs {e2}");
}

#[test]
fn direct_reciprocal_forces_match_numeric_gradient() {
    let sys = neutral_system(55);
    let alpha = 3.8;
    let h = 1e-5;
    let (_, f) = direct_total(&sys, alpha);
    for atom in [0usize, 5, 11] {
        let energy_at = |dx: f64| {
            let mut s =
                Neutral { x: sys.x.clone(), y: sys.y.clone(), z: sys.z.clone(), q: sys.q.clone() };
            s.x[atom] += dx;
            direct_total(&s, alpha).0
        };
        let numeric = -(energy_at(h) - energy_at(-h)) / (2.0 * h);
        assert!(
            (f.0[atom] - numeric).abs() < 1e-5,
            "atom {atom}: analytic {} numeric {numeric}",
            f.0[atom]
        );
    }
}

#[test]
fn pme_reciprocal_energy_and_forces_match_direct_ewald() {
    let sys = neutral_system(1234);
    let alpha = 3.8;
    let pos = Positions::new(&sys.x, &sys.y, &sys.z);
    let box_len = [EDGE; 3];

    let mut fd = zeros();
    let e_direct = direct_reciprocal(
        pos,
        &sys.q,
        box_len,
        alpha,
        10,
        &mut ForceBuffers::new(&mut fd.0, &mut fd.1, &mut fd.2),
    );

    let params = PmeParams { alpha, grid: [32, 32, 32], order: 6 };
    let pme = Pme::new(box_len, params);
    let mut fp = zeros();
    let e_pme =
        pme.reciprocal(pos, &sys.q, &mut ForceBuffers::new(&mut fp.0, &mut fp.1, &mut fp.2));

    let e_err = rel_err(e_pme, e_direct);
    let f_err = rms_rel(&fp, &fd);
    println!("PME reciprocal: E_direct={e_direct:.8} E_pme={e_pme:.8} rel_err={e_err:.2e} rms_force_rel_err={f_err:.2e}");
    assert!(e_err < 1e-5, "reciprocal energy relative error {e_err:.2e}");
    assert!(f_err < 1e-4, "reciprocal force RMS relative error {f_err:.2e}");
}

#[test]
fn pme_error_decreases_with_order_and_grid() {
    let sys = neutral_system(1234);
    let alpha = 3.8;
    let pos = Positions::new(&sys.x, &sys.y, &sys.z);
    let box_len = [EDGE; 3];
    let mut fd = zeros();
    let e_direct = direct_reciprocal(
        pos,
        &sys.q,
        box_len,
        alpha,
        10,
        &mut ForceBuffers::new(&mut fd.0, &mut fd.1, &mut fd.2),
    );
    let err = |grid: usize, order: usize| {
        let pme = Pme::new(box_len, PmeParams { alpha, grid: [grid; 3], order });
        let mut fp = zeros();
        let e =
            pme.reciprocal(pos, &sys.q, &mut ForceBuffers::new(&mut fp.0, &mut fp.1, &mut fp.2));
        rel_err(e, e_direct)
    };
    let coarse = err(16, 4);
    let fine = err(32, 6);
    println!("PME energy error: grid16/order4={coarse:.2e} grid32/order6={fine:.2e}");
    assert!(fine < coarse, "finer grid/higher order should be more accurate");
}

#[test]
fn full_pme_pipeline_matches_direct_ewald() {
    // Real-space (min image, cutoff L/2) + PME reciprocal + self energy, against
    // the converged direct Ewald total. Real-space truncation at L/2 leaves
    // ~erfc(3.8) per pair, so the energy tolerance is set accordingly.
    let sys = neutral_system(4321);
    let alpha = 3.8;
    let (e_direct, fd) = direct_total(&sys, alpha);

    let pos = Positions::new(&sys.x, &sys.y, &sys.z);
    let params = AtomParams { sigma: &[0.0; N], epsilon: &[0.0; N], charge: &sys.q };
    let cfg = NonbondedConfig { cutoff: 1.0, alpha, box_len: Some([EDGE; 3]) };
    let pairs: Vec<(u32, u32)> =
        (0..N as u32).flat_map(|i| (i + 1..N as u32).map(move |j| (i, j))).collect();
    let mut fr = zeros();
    let e_real = pair_list_scalar(
        pos,
        params,
        &pairs,
        &cfg,
        &mut ForceBuffers::new(&mut fr.0, &mut fr.1, &mut fr.2),
    );

    let pme = Pme::new([EDGE; 3], PmeParams { alpha, grid: [32; 3], order: 6 });
    let mut fk = zeros();
    let e_rec =
        pme.reciprocal(pos, &sys.q, &mut ForceBuffers::new(&mut fk.0, &mut fk.1, &mut fk.2));
    let e_self = pme.self_energy(&sys.q);
    let e_total = e_real.coulomb + e_rec + e_self;

    let mut ftot = fr;
    for k in 0..N {
        ftot.0[k] += fk.0[k];
        ftot.1[k] += fk.1[k];
        ftot.2[k] += fk.2[k];
    }
    let e_err = (e_total - e_direct).abs();
    let f_err = rms_rel(&ftot, &fd);
    println!("full PME: E_direct={e_direct:.8} E_pme={e_total:.8} abs_err={e_err:.2e} rms_force_rel_err={f_err:.2e}");
    assert!(e_err < 1e-4, "total energy absolute error {e_err:.2e}");
    assert!(f_err < 1e-4, "total force RMS relative error {f_err:.2e}");

    // The SIMD and cluster real-space kernels must agree with the scalar one here too.
    let mut fs = zeros();
    let e_simd = pair_list_simd(
        pos,
        params,
        &pairs,
        &cfg,
        &mut ForceBuffers::new(&mut fs.0, &mut fs.1, &mut fs.2),
    );
    assert!((e_simd.coulomb - e_real.coulomb).abs() < 1e-12);
    let set = ClusterSet::from_atoms(pos, params);
    let cp = build_cluster_pairs(&set, cfg.cutoff, cfg.box_len);
    let mut cf = set.zero_forces();
    let e_cl = cluster_pair_kernel(&set, &cp, &cfg, &mut cf);
    assert!((e_cl.coulomb - e_real.coulomb).abs() < 1e-12);
}

#[test]
fn parameter_selection_hits_accuracy_target() {
    let tol = 1e-5;
    let rc = 1.0;
    let alpha = alpha_for_tolerance(rc, tol);
    assert!((libm::erfc(alpha * rc) - tol).abs() < 1e-10 * tol.max(1.0));

    let p = PmeParams::from_accuracy([EDGE; 3], rc, tol, 6);
    assert!((p.alpha - alpha).abs() < 1e-12);
    assert_eq!(p.order, 6);
    for n in p.grid {
        assert!(n.is_power_of_two());
        assert!(EDGE / n as f64 <= 0.12 + 1e-12, "spacing too coarse: {}", EDGE / n as f64);
    }
}

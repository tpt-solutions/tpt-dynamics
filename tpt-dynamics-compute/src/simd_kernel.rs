//! SIMD pair-list kernel built on `tpt-simd` 4-lane `f64` vectors.
//!
//! Four pairs are processed per iteration: separations, minimum-image wrapping,
//! cutoff masks and the Lennard-Jones polynomial run in vector registers. The
//! real-space Ewald `erfc` and `exp` are evaluated per lane with `libm` (there is
//! no vector `erfc` in `tpt-simd`), and force scatter into the atom arrays is
//! scalar because pairs share atoms.

use tpt_simd::{F64x4, SimdMask};

use crate::pairwise::{coulomb_ewald_pair, lorentz_berthelot};
use crate::types::{AtomParams, Energies, ForceBuffers, NonbondedConfig, Positions};

/// Number of lanes processed together.
pub const LANES: usize = 4;

/// Minimum-image wrapping applied lane-wise to a vector of separations.
#[inline]
pub(crate) fn min_image_v(d: F64x4, edge: f64) -> F64x4 {
    let l = F64x4::splat(edge);
    d - l * (d / l).round()
}

/// Sums the lanes of a vector.
#[inline]
pub(crate) fn hsum(v: F64x4) -> f64 {
    let a: [f64; LANES] = v.into();
    a.iter().sum()
}

/// Evaluates LJ and real-space Ewald terms for up to four pairs at once.
///
/// Inputs are per-lane: squared separation `r2`, a validity `mask` (lanes that
/// must contribute), combined `sigma`/`epsilon` and the charge product `qq`.
/// Returns `(e_lj, e_coulomb, f_over_r)` with unmasked lanes set to zero.
#[inline]
pub(crate) fn lane_interaction(
    r2: F64x4,
    mask: SimdMask<f64, LANES>,
    sigma: F64x4,
    epsilon: F64x4,
    qq: F64x4,
    alpha: f64,
) -> (F64x4, F64x4, F64x4) {
    let zero = F64x4::splat(0.0);
    // Masked lanes get r2 = 1 so that no division produces inf/NaN.
    let r2s = mask.select(r2, F64x4::splat(1.0));

    let sr2 = sigma * sigma / r2s;
    let sr6 = sr2 * sr2 * sr2;
    let sr12 = sr6 * sr6;
    let e_lj = mask.select(F64x4::splat(4.0) * epsilon * (sr12 - sr6), zero);
    let f_lj = F64x4::splat(24.0) * epsilon * (F64x4::splat(2.0) * sr12 - sr6) / r2s;

    let r2a: [f64; LANES] = r2s.into();
    let qqa: [f64; LANES] = qq.into();
    let ma = mask.to_array();
    let mut ec = [0.0; LANES];
    let mut fc = [0.0; LANES];
    for lane in 0..LANES {
        if ma[lane] {
            let (e, f) = coulomb_ewald_pair(r2a[lane], qqa[lane], alpha);
            ec[lane] = e;
            fc[lane] = f;
        }
    }
    let e_c = F64x4::from(ec);
    let f_total = mask.select(f_lj + F64x4::from(fc), zero);
    (e_lj, e_c, f_total)
}

/// SIMD version of [`crate::scalar::pair_list_scalar`].
///
/// Same inputs, same outputs (forces are added, energies returned). Results agree
/// with the scalar kernel to round-off; see the crate tests.
pub fn pair_list_simd(
    pos: Positions<'_>,
    params: AtomParams<'_>,
    pairs: &[(u32, u32)],
    cfg: &NonbondedConfig,
    forces: &mut ForceBuffers<'_>,
) -> Energies {
    let rc2 = cfg.cutoff * cfg.cutoff;
    let rc2v = F64x4::splat(rc2);
    let zero = F64x4::splat(0.0);
    let mut en = Energies::default();

    for chunk in pairs.chunks(LANES) {
        let mut dx = [0.0; LANES];
        let mut dy = [0.0; LANES];
        let mut dz = [0.0; LANES];
        let mut sig = [0.0; LANES];
        let mut eps = [0.0; LANES];
        let mut qq = [0.0; LANES];
        let mut valid = [false; LANES];
        for (lane, &(i, j)) in chunk.iter().enumerate() {
            let (iu, ju) = (i as usize, j as usize);
            // Raw separations; minimum image is applied to the whole vector below.
            dx[lane] = pos.x[iu] - pos.x[ju];
            dy[lane] = pos.y[iu] - pos.y[ju];
            dz[lane] = pos.z[iu] - pos.z[ju];
            let (s, e) = lorentz_berthelot(
                params.sigma[iu],
                params.epsilon[iu],
                params.sigma[ju],
                params.epsilon[ju],
            );
            sig[lane] = s;
            eps[lane] = e;
            qq[lane] = params.charge[iu] * params.charge[ju];
            valid[lane] = true;
        }

        let (mut dxv, mut dyv, mut dzv) = (F64x4::from(dx), F64x4::from(dy), F64x4::from(dz));
        if let Some(l) = cfg.box_len {
            dxv = min_image_v(dxv, l[0]);
            dyv = min_image_v(dyv, l[1]);
            dzv = min_image_v(dzv, l[2]);
        }
        let r2 = dxv * dxv + dyv * dyv + dzv * dzv;
        let mask = r2.simd_lt(rc2v) & r2.simd_gt(zero) & SimdMask::<f64, LANES>::from_array(valid);

        let (e_lj, e_c, f) = lane_interaction(
            r2,
            mask,
            F64x4::from(sig),
            F64x4::from(eps),
            F64x4::from(qq),
            cfg.alpha,
        );
        en.lj += hsum(e_lj);
        en.coulomb += hsum(e_c);

        let fxa: [f64; LANES] = (f * dxv).into();
        let fya: [f64; LANES] = (f * dyv).into();
        let fza: [f64; LANES] = (f * dzv).into();
        let ma = mask.to_array();
        for (lane, &(i, j)) in chunk.iter().enumerate() {
            if ma[lane] {
                let (iu, ju) = (i as usize, j as usize);
                forces.fx[iu] += fxa[lane];
                forces.fy[iu] += fya[lane];
                forces.fz[iu] += fza[lane];
                forces.fx[ju] -= fxa[lane];
                forces.fy[ju] -= fya[lane];
                forces.fz[ju] -= fza[lane];
            }
        }
    }
    en
}

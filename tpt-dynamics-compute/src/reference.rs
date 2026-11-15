//! Brute-force Ewald summation used as the accuracy reference for PME.
//!
//! These routines are O(N * K^3) and O(N^2 * images) and are meant for validation
//! on small systems, not for production use.
//!
//! * [`direct_reciprocal`] sums the reciprocal-space Ewald series over integer
//!   wave-vectors `n` with `|n_axis| <= kmax` (`k = 2 pi n / L`).
//! * [`direct_real_images`] sums the real-space `erfc` term over periodic images
//!   `|n_axis| <= nimg`.
//!
//! Energies use the Coulomb constant 1 and include the self-energy term via
//! [`direct_self_energy`].

use crate::pairwise::coulomb_ewald_pair;
use crate::types::{ForceBuffers, Positions};

/// Self-energy `-alpha/sqrt(pi) * sum q^2`.
pub fn direct_self_energy(charges: &[f64], alpha: f64) -> f64 {
    -alpha / std::f64::consts::PI.sqrt() * charges.iter().map(|q| q * q).sum::<f64>()
}

/// Reciprocal-space Ewald energy and forces by direct summation over `|n| <= kmax`.
///
/// `E = (2 pi / V) sum_{k != 0} exp(-k^2 / 4 alpha^2) / k^2 |S(k)|^2`, with
/// `S(k) = sum_j q_j exp(i k . r_j)`. Forces are added to `forces`.
pub fn direct_reciprocal(
    pos: Positions<'_>,
    charges: &[f64],
    box_len: [f64; 3],
    alpha: f64,
    kmax: i32,
    forces: &mut ForceBuffers<'_>,
) -> f64 {
    let pi = std::f64::consts::PI;
    let volume = box_len[0] * box_len[1] * box_len[2];
    let n = pos.len();
    let mut energy = 0.0;
    for nx in -kmax..=kmax {
        for ny in -kmax..=kmax {
            for nz in -kmax..=kmax {
                if nx == 0 && ny == 0 && nz == 0 {
                    continue;
                }
                let k = [
                    2.0 * pi * f64::from(nx) / box_len[0],
                    2.0 * pi * f64::from(ny) / box_len[1],
                    2.0 * pi * f64::from(nz) / box_len[2],
                ];
                let k2 = k[0] * k[0] + k[1] * k[1] + k[2] * k[2];
                let g = (-k2 / (4.0 * alpha * alpha)).exp() / k2;
                // Phases e^{i k.r_j}.
                let mut s = (0.0, 0.0);
                let mut phase = vec![(0.0, 0.0); n];
                for j in 0..n {
                    let ph = k[0] * pos.x[j] + k[1] * pos.y[j] + k[2] * pos.z[j];
                    let (sn, cs) = (ph.sin(), ph.cos());
                    phase[j] = (cs, sn);
                    s.0 += charges[j] * cs;
                    s.1 += charges[j] * sn;
                }
                energy += (2.0 * pi / volume) * g * (s.0 * s.0 + s.1 * s.1);
                // F_i = (4 pi q_i / V) * g * k * Im(S^* e^{i k.r_i})
                for i in 0..n {
                    let (cs, sn) = phase[i];
                    // S^* e^{i k r_i} = (s0 - i s1)(cs + i sn): imaginary part = s0 sn - s1 cs.
                    let im = s.0 * sn - s.1 * cs;
                    let f = 4.0 * pi / volume * charges[i] * g * im;
                    forces.fx[i] += f * k[0];
                    forces.fy[i] += f * k[1];
                    forces.fz[i] += f * k[2];
                }
            }
        }
    }
    energy
}

/// Real-space Ewald energy over periodic images `|n| <= nimg`, with forces.
///
/// Sums `qq erfc(alpha r)/r` over all pairs and images. Self-images contribute energy
/// only (their forces vanish). Pairs use the full image sum, so there is no cutoff.
pub fn direct_real_images(
    pos: Positions<'_>,
    charges: &[f64],
    box_len: [f64; 3],
    alpha: f64,
    nimg: i32,
    forces: &mut ForceBuffers<'_>,
) -> f64 {
    let n = pos.len();
    let mut energy = 0.0;
    for i in 0..n {
        for j in 0..n {
            for nx in -nimg..=nimg {
                for ny in -nimg..=nimg {
                    for nz in -nimg..=nimg {
                        // Unordered pairs for i != j; ordered (with the 1/2) for self-images.
                        if i == j && nx == 0 && ny == 0 && nz == 0 {
                            continue;
                        }
                        if i > j {
                            continue;
                        }
                        let d = [
                            pos.x[i] - pos.x[j] + f64::from(nx) * box_len[0],
                            pos.y[i] - pos.y[j] + f64::from(ny) * box_len[1],
                            pos.z[i] - pos.z[j] + f64::from(nz) * box_len[2],
                        ];
                        let r2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
                        let qq = charges[i] * charges[j];
                        if i == j {
                            // Self-image: half of each ordered (n, -n) pair; no force.
                            let (e, _) = coulomb_ewald_pair(r2, qq, alpha);
                            energy += 0.5 * e;
                        } else {
                            let (e, f_over_r) = coulomb_ewald_pair(r2, qq, alpha);
                            energy += e;
                            forces.fx[i] += f_over_r * d[0];
                            forces.fy[i] += f_over_r * d[1];
                            forces.fz[i] += f_over_r * d[2];
                            forces.fx[j] -= f_over_r * d[0];
                            forces.fy[j] -= f_over_r * d[1];
                            forces.fz[j] -= f_over_r * d[2];
                        }
                    }
                }
            }
        }
    }
    energy
}

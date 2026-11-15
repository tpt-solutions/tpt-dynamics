//! Energy minimisation: adaptive steepest descent and limited-memory BFGS.
//!
//! Both minimisers move atomic positions to reduce the potential energy from
//! the force field in a [`ForceEvaluator`]. Positions are not wrapped during
//! minimisation, so displacement vectors stay continuous. The minimum-image
//! convention is used for all pair and bonded terms, and the final positions are
//! wrapped into the box.
//!
//! * [`steepest_descent`]: moves the atom with the largest force by an adaptive
//!   step `h`. A step is accepted only if the energy decreases; then `h` grows
//!   by 1.2. Otherwise it is rejected and `h` is multiplied by 0.2. Accepted
//!   energies therefore decrease monotonically.
//! * [`lbfgs`]: limited-memory BFGS with the standard two-loop recursion,
//!   backtracking Armijo line search, and a cap on the largest atomic
//!   displacement per step.
//!
//! Convergence means the largest force component is below
//! [`MinimizeSettings::force_tolerance`].

use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::integrator::ForceEvaluator;
use crate::system::System;

/// Minimiser settings.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MinimizeSettings {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence threshold on the largest force component (kJ mol⁻¹ nm⁻¹).
    pub force_tolerance: f64,
    /// Initial steepest-descent step length (nm).
    pub initial_step: f64,
    /// Largest atomic displacement allowed in one step (nm).
    pub max_displacement: f64,
    /// Number of correction pairs kept by L-BFGS.
    pub memory: usize,
}

impl Default for MinimizeSettings {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            force_tolerance: 10.0,
            initial_step: 0.01,
            max_displacement: 0.1,
            memory: 5,
        }
    }
}

/// Outcome of a minimisation.
#[derive(Clone, Debug, PartialEq)]
pub struct MinimizeReport {
    /// Iterations performed.
    pub iterations: usize,
    /// Potential energy before minimisation.
    pub initial_energy: f64,
    /// Potential energy after minimisation.
    pub final_energy: f64,
    /// Largest force component after minimisation.
    pub max_force: f64,
    /// True if `max_force` is below the tolerance.
    pub converged: bool,
    /// Energy after each accepted iteration (starting with the initial energy).
    pub history: Vec<f64>,
}

fn validate(settings: &MinimizeSettings) -> Result<(), SimError> {
    if !(settings.force_tolerance > 0.0
        && settings.initial_step > 0.0
        && settings.max_displacement > 0.0)
    {
        return Err(SimError::Config("minimiser tolerances and steps must be > 0".into()));
    }
    if settings.memory == 0 {
        return Err(SimError::Config("L-BFGS memory must be >= 1".into()));
    }
    Ok(())
}

fn max_force(forces: &[Vec3]) -> f64 {
    forces.iter().flat_map(|f| [f.x.abs(), f.y.abs(), f.z.abs()]).fold(0.0, f64::max)
}

fn finish(system: &mut System, eval: &mut ForceEvaluator) -> Result<f64, SimError> {
    let sb = system.sim_box;
    for p in &mut system.positions {
        *p = sb.wrap(*p);
    }
    Ok(eval.evaluate(system)?.potential)
}

/// Adaptive steepest-descent minimisation.
///
/// # Errors
/// Returns [`SimError::Config`] for invalid settings and propagates force errors.
pub fn steepest_descent(
    system: &mut System,
    eval: &mut ForceEvaluator,
    settings: &MinimizeSettings,
) -> Result<MinimizeReport, SimError> {
    validate(settings)?;
    let out = eval.evaluate(system)?;
    let mut energy = out.potential;
    let initial_energy = energy;
    let mut history = vec![energy];
    let mut step = settings.initial_step;
    let mut iterations = 0;
    let mut converged = false;
    while iterations < settings.max_iterations {
        let fmax = max_force(&system.forces);
        if fmax < settings.force_tolerance {
            converged = true;
            break;
        }
        if step < 1e-14 {
            break;
        }
        iterations += 1;
        // Largest atom moves by `step`, the others proportionally to their force.
        let scale = step.min(settings.max_displacement) / fmax;
        let old_positions = system.positions.clone();
        let old_forces = system.forces.clone();
        for (p, f) in system.positions.iter_mut().zip(&old_forces) {
            *p += *f * scale;
        }
        let trial = eval.evaluate(system)?.potential;
        if trial < energy {
            energy = trial;
            step *= 1.2;
            history.push(energy);
        } else {
            system.positions = old_positions;
            system.forces = old_forces;
            step *= 0.2;
            // Restore the list and energies for the accepted configuration.
            eval.evaluate(system)?;
        }
    }
    let final_energy = finish(system, eval)?;
    Ok(MinimizeReport {
        iterations,
        initial_energy,
        final_energy,
        max_force: max_force(&system.forces),
        converged,
        history,
    })
}

fn flatten(v: &[Vec3]) -> Vec<f64> {
    v.iter().flat_map(|p| [p.x, p.y, p.z]).collect()
}

fn unflatten(x: &[f64]) -> Vec<Vec3> {
    x.as_chunks::<3>().0.iter().map(|c| Vec3::new(c[0], c[1], c[2])).collect()
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

/// L-BFGS minimisation.
///
/// # Errors
/// Returns [`SimError::Config`] for invalid settings and propagates force errors.
pub fn lbfgs(
    system: &mut System,
    eval: &mut ForceEvaluator,
    settings: &MinimizeSettings,
) -> Result<MinimizeReport, SimError> {
    validate(settings)?;
    let out = eval.evaluate(system)?;
    let mut energy = out.potential;
    let initial_energy = energy;
    let mut history = vec![energy];
    let m = settings.memory;
    let mut s_hist: Vec<Vec<f64>> = Vec::new();
    let mut y_hist: Vec<Vec<f64>> = Vec::new();
    let mut x = flatten(&system.positions);
    let mut g: Vec<f64> = flatten(&system.forces).iter().map(|f| -f).collect();
    let mut iterations = 0;
    let mut converged = false;
    while iterations < settings.max_iterations {
        if max_force(&system.forces) < settings.force_tolerance {
            converged = true;
            break;
        }
        iterations += 1;
        // Two-loop recursion for d = −H g.
        let mut q = g.clone();
        let k = s_hist.len();
        let mut alpha = vec![0.0; k];
        for i in (0..k).rev() {
            let rho = 1.0 / dot(&y_hist[i], &s_hist[i]);
            alpha[i] = rho * dot(&s_hist[i], &q);
            for (qj, yj) in q.iter_mut().zip(&y_hist[i]) {
                *qj -= alpha[i] * yj;
            }
        }
        let gamma = if k > 0 {
            dot(&s_hist[k - 1], &y_hist[k - 1]) / dot(&y_hist[k - 1], &y_hist[k - 1])
        } else {
            // First step: a conservative length based on the step settings.
            let gn = dot(&g, &g).sqrt().max(1e-300);
            settings.initial_step / gn
        };
        for qj in &mut q {
            *qj *= gamma;
        }
        for i in 0..k {
            let rho = 1.0 / dot(&y_hist[i], &s_hist[i]);
            let beta = rho * dot(&y_hist[i], &q);
            for (qj, sj) in q.iter_mut().zip(&s_hist[i]) {
                *qj += (alpha[i] - beta) * sj;
            }
        }
        let mut d: Vec<f64> = q.iter().map(|v| -v).collect();
        // Not a descent direction (possible with rounding): fall back to steepest descent.
        if dot(&d, &g) >= 0.0 {
            d = g
                .iter()
                .map(|v| -v * settings.initial_step / dot(&g, &g).sqrt().max(1e-300))
                .collect();
            s_hist.clear();
            y_hist.clear();
        }
        // Cap the largest atomic displacement.
        let dmax = d.iter().fold(0.0_f64, |a, v| a.max(v.abs()));
        let mut step =
            if dmax > settings.max_displacement { settings.max_displacement / dmax } else { 1.0 };
        // Backtracking Armijo line search.
        let slope = dot(&g, &d);
        let mut accepted = false;
        for _ in 0..40 {
            let trial: Vec<f64> = x.iter().zip(&d).map(|(xi, di)| xi + step * di).collect();
            system.positions = unflatten(&trial);
            let e = eval.evaluate(system)?.potential;
            if e <= energy + 1e-4 * step * slope {
                accepted = true;
                let g_new: Vec<f64> = flatten(&system.forces).iter().map(|f| -f).collect();
                let s: Vec<f64> = trial.iter().zip(&x).map(|(a, b)| a - b).collect();
                let y: Vec<f64> = g_new.iter().zip(&g).map(|(a, b)| a - b).collect();
                if dot(&s, &y) > 1e-14 {
                    if s_hist.len() == m {
                        s_hist.remove(0);
                        y_hist.remove(0);
                    }
                    s_hist.push(s);
                    y_hist.push(y);
                }
                x = trial;
                g = g_new;
                energy = e;
                history.push(energy);
                break;
            }
            step *= 0.5;
        }
        if !accepted {
            // No decrease possible along this direction: restore and stop.
            system.positions = unflatten(&x);
            eval.evaluate(system)?;
            break;
        }
    }
    let final_energy = finish(system, eval)?;
    Ok(MinimizeReport {
        iterations,
        initial_energy,
        final_energy,
        max_force: max_force(&system.forces),
        converged,
        history,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::forcefield::{LennardJones, Modifier};
    use crate::pbc::SimBox;
    use crate::system::NewAtom;

    fn lj_system(positions: &[Vec3], l: f64) -> System {
        let mut sys = System::new(SimBox::new(l, l, l).unwrap());
        for &p in positions {
            sys.push_atom(NewAtom {
                name: "AR".into(),
                residue_name: "LJ".into(),
                residue_id: 1,
                element: "Ar".into(),
                position: p + Vec3::new(l / 2.0, l / 2.0, l / 2.0),
                velocity: Vec3::ZERO,
                mass: 1.0,
                charge: 0.0,
                atom_type: 0,
            });
        }
        sys
    }

    fn evaluator(cutoff: f64) -> ForceEvaluator {
        // No shift: the cluster energies then match the textbook LJ values.
        let ff = LennardJones::single_type(1.0, 1.0, cutoff, Modifier::None).unwrap();
        ForceEvaluator::new(Box::new(ff), 0.3).unwrap()
    }

    fn cluster13() -> Vec<Vec3> {
        // Icosahedron-like start: a centre atom with a shell of 12 at distance ~1.1.
        let mut pts = vec![Vec3::ZERO];
        let phi = f64::midpoint(1.0, 5f64.sqrt());
        let raw = [
            (0.0, 1.0, phi),
            (0.0, -1.0, phi),
            (0.0, 1.0, -phi),
            (0.0, -1.0, -phi),
            (1.0, phi, 0.0),
            (-1.0, phi, 0.0),
            (1.0, -phi, 0.0),
            (-1.0, -phi, 0.0),
            (phi, 0.0, 1.0),
            (-phi, 0.0, 1.0),
            (phi, 0.0, -1.0),
            (-phi, 0.0, -1.0),
        ];
        let norm = (1.0 + phi * phi).sqrt();
        for (x, y, z) in raw {
            pts.push(Vec3::new(x, y, z) * (1.1 / norm));
        }
        pts
    }

    #[test]
    fn steepest_descent_energy_decreases_monotonically() {
        let mut sys = lj_system(&cluster13(), 40.0);
        let mut eval = evaluator(5.0);
        let settings =
            MinimizeSettings { max_iterations: 400, force_tolerance: 1e-3, ..Default::default() };
        let rep = steepest_descent(&mut sys, &mut eval, &settings).unwrap();
        assert!(rep.history.windows(2).all(|w| w[1] < w[0]), "not monotone");
        assert!(rep.final_energy < rep.initial_energy);
        assert!((rep.final_energy - eval.evaluate(&mut sys).unwrap().potential).abs() < 1e-12);
    }

    #[test]
    fn lj_dimer_converges_to_two_to_one_sixth_sigma() {
        let r_min = 2f64.powf(1.0 / 6.0);
        for use_lbfgs in [false, true] {
            let mut sys = lj_system(&[Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)], 20.0);
            let mut eval = evaluator(5.0);
            // Steepest descent is slow near the minimum, hence the generous limits.
            let settings = MinimizeSettings {
                force_tolerance: 1e-4,
                max_iterations: 20_000,
                ..Default::default()
            };
            let rep = if use_lbfgs {
                lbfgs(&mut sys, &mut eval, &settings).unwrap()
            } else {
                steepest_descent(&mut sys, &mut eval, &settings).unwrap()
            };
            assert!(rep.converged, "lbfgs={use_lbfgs} did not converge: {rep:?}");
            let r = (sys.positions[1] - sys.positions[0]).norm();
            assert!((r - r_min).abs() < 1e-3, "lbfgs={use_lbfgs}: r = {r}");
            assert!((rep.final_energy + 1.0).abs() < 1e-6, "E = {}", rep.final_energy);
        }
    }

    #[test]
    fn lbfgs_finds_thirteen_atom_lj_global_minimum() {
        // Global minimum of LJ13 (ε = σ = 1, no cutoff): E = −44.326801.
        let mut sys = lj_system(&cluster13(), 40.0);
        let mut eval = evaluator(5.0);
        let settings =
            MinimizeSettings { force_tolerance: 1e-5, max_iterations: 5000, ..Default::default() };
        let rep = lbfgs(&mut sys, &mut eval, &settings).unwrap();
        assert!(rep.converged, "max force {}", rep.max_force);
        assert!((rep.final_energy + 44.326_801).abs() < 1e-4, "E = {}", rep.final_energy);
    }

    #[test]
    fn invalid_settings_are_rejected() {
        let mut sys = lj_system(&[Vec3::ZERO, Vec3::new(1.0, 0.0, 0.0)], 20.0);
        let mut eval = evaluator(5.0);
        let bad = MinimizeSettings { memory: 0, ..Default::default() };
        assert!(lbfgs(&mut sys, &mut eval, &bad).is_err());
    }
}

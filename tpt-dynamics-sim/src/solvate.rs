//! Solvation and ion placement.
//!
//! * [`place_water_grid`] fills a box with placeholder three-site water molecules
//!   on a regular grid. The geometry is the TIP3P-like O–H and H–H distances,
//!   but the charges are zero and the molecules have no force-field
//!   parameters attached. The output is for geometry tests and as input for
//!   the ion placement below; it is not a equilibrated solvent.
//! * [`replace_waters_with_ions`] removes randomly chosen waters and puts Na⁺ and
//!   Cl⁻ ions at their oxygen sites. Choices are deterministic for a given seed,
//!   and a candidate is rejected if it lies closer than `min_separation` to an
//!   ion already placed. The number of Na⁺ and Cl⁻ ions is equal, so the net
//!   charge is unchanged.

use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::pbc::SimBox;
use crate::rng::Rng;
use crate::system::{NewAtom, System};

/// O–H distance (nm) of the placeholder water geometry.
pub const WATER_OH: f64 = 0.09572;
/// H–H distance (nm) of the placeholder water geometry.
pub const WATER_HH: f64 = 0.15139;
/// Mass of oxygen (amu).
pub const MASS_O: f64 = 15.9994;
/// Mass of hydrogen (amu).
pub const MASS_H: f64 = 1.008;

/// Fills `sim_box` with placeholder water on a grid of cells with edge at most
/// `spacing` (nm). Each cell holds one molecule whose oxygen is at the cell
/// centre; molecules lie in the xy plane, bisector along +y.
///
/// # Errors
/// Returns [`SimError::Invalid`] if `spacing` is not positive or is so large that
/// no molecule fits.
pub fn place_water_grid(sim_box: SimBox, spacing: f64) -> Result<System, SimError> {
    place_water_grid_geometry(sim_box, spacing, WATER_OH, WATER_HH)
}

/// As [`place_water_grid`], with an explicit O–H and H–H distance (nm), for
/// example the SPC/E geometry from [`crate::water::WaterModel::SpcE`].
///
/// # Errors
/// Returns [`SimError::Invalid`] if `spacing`, `d_oh` or `d_hh` is not positive
/// and finite, or if `d_hh` exceeds `2 d_oh`.
pub fn place_water_grid_geometry(
    sim_box: SimBox,
    spacing: f64,
    d_oh: f64,
    d_hh: f64,
) -> Result<System, SimError> {
    if !(spacing.is_finite() && spacing > 0.0) {
        return Err(SimError::Invalid("grid spacing must be > 0".into()));
    }
    if !(d_oh.is_finite() && d_oh > 0.0 && d_hh.is_finite() && d_hh > 0.0 && d_hh < 2.0 * d_oh) {
        return Err(SimError::Invalid("water O-H and H-H distances must be positive".into()));
    }
    let l = sim_box.lengths;
    let n = [
        ((l.x / spacing).floor() as usize).max(1),
        ((l.y / spacing).floor() as usize).max(1),
        ((l.z / spacing).floor() as usize).max(1),
    ];
    let mut sys = System::new(sim_box);
    let half = 0.5 * d_hh;
    let height = (d_oh * d_oh - half * half).sqrt();
    let mut resid = 0;
    for ix in 0..n[0] {
        for iy in 0..n[1] {
            for iz in 0..n[2] {
                resid += 1;
                let centre = Vec3::new(
                    (ix as f64 + 0.5) * l.x / n[0] as f64,
                    (iy as f64 + 0.5) * l.y / n[1] as f64,
                    (iz as f64 + 0.5) * l.z / n[2] as f64,
                );
                let geometry = [
                    (centre, "OW", "O", MASS_O),
                    (centre + Vec3::new(half, -height, 0.0), "HW1", "H", MASS_H),
                    (centre + Vec3::new(-half, -height, 0.0), "HW2", "H", MASS_H),
                ];
                for (pos, name, element, mass) in geometry {
                    sys.push_atom(NewAtom {
                        name: name.into(),
                        residue_name: "SOL".into(),
                        residue_id: resid,
                        element: element.into(),
                        position: pos,
                        velocity: Vec3::ZERO,
                        mass,
                        charge: 0.0,
                        atom_type: 0,
                    });
                }
            }
        }
    }
    Ok(sys)
}

/// Summary of an ion placement.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct IonPlacement {
    /// Na⁺ ions placed (equal to the number of Cl⁻ ions).
    pub sodium: usize,
    /// Cl⁻ ions placed.
    pub chloride: usize,
    /// Waters removed (three atoms each).
    pub waters_removed: usize,
}

/// Replaces `n_pairs` waters with Na⁺ and `n_pairs` waters with Cl⁻.
///
/// The system must consist of three-atom water blocks (as made by
/// [`place_water_grid`]): each block is an oxygen followed by two hydrogens.
/// Candidate oxygens are visited in an order drawn from a seeded generator, and
/// a candidate is skipped if it lies within `min_separation` (nm) of an ion
/// already placed.
///
/// # Errors
/// Returns [`SimError::Invalid`] if the system is not a whole number of water
/// blocks, or if too few candidates satisfy the separation rule.
pub fn replace_waters_with_ions(
    system: &mut System,
    n_pairs: usize,
    seed: u64,
    min_separation: f64,
) -> Result<IonPlacement, SimError> {
    if !system.len().is_multiple_of(3) {
        return Err(SimError::Invalid("system is not a whole number of water blocks".into()));
    }
    if !(min_separation.is_finite() && min_separation >= 0.0) {
        return Err(SimError::Invalid("minimum ion separation must be >= 0".into()));
    }
    let n_water = system.len() / 3;
    let sim_box = system.sim_box;
    let mut order: Vec<usize> = (0..n_water).collect();
    let mut rng = Rng::seed_from_u64(seed);
    // Fisher-Yates shuffle, deterministic for the seed.
    for i in (1..order.len()).rev() {
        let j = (rng.next_u64() % (i as u64 + 1)) as usize;
        order.swap(i, j);
    }

    let mut chosen: Vec<usize> = Vec::new();
    let mut sites: Vec<Vec3> = Vec::new();
    let want = 2 * n_pairs;
    for &w in &order {
        if chosen.len() == want {
            break;
        }
        let o = system.positions[3 * w];
        if sites.iter().all(|&s| sim_box.min_image(o - s).norm() >= min_separation) {
            chosen.push(w);
            sites.push(o);
        }
    }
    if chosen.len() < want {
        return Err(SimError::Invalid(format!(
            "only {} of {want} ion sites satisfy the separation rule",
            chosen.len()
        )));
    }

    // Rebuild the atom list: keep waters not chosen, add ions at the chosen oxygens.
    let mut keep = vec![true; system.len()];
    for &w in &chosen {
        keep[3 * w] = false;
        keep[3 * w + 1] = false;
        keep[3 * w + 2] = false;
    }
    let ion_sites: Vec<(usize, Vec3)> = chosen.iter().copied().zip(sites.iter().copied()).collect();
    let mut ions = Vec::new();
    for (k, (w, pos)) in ion_sites.into_iter().enumerate() {
        let (name, element, mass, charge) =
            if k < n_pairs { ("NA", "Na", 22.990, 1.0) } else { ("CL", "Cl", 35.45, -1.0) };
        ions.push((
            w,
            NewAtom {
                name: name.into(),
                residue_name: name.into(),
                residue_id: 0,
                element: element.into(),
                position: pos,
                velocity: Vec3::ZERO,
                mass,
                charge,
                atom_type: 0,
            },
        ));
    }
    system.retain(&keep);
    for (_, atom) in ions {
        system.push_atom(atom);
    }
    Ok(IonPlacement { sodium: n_pairs, chloride: n_pairs, waters_removed: want })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn box10() -> SimBox {
        SimBox::new(3.0, 3.0, 3.0).unwrap()
    }

    #[test]
    fn grid_places_one_molecule_per_cell_with_correct_geometry() {
        let sys = place_water_grid(box10(), 1.0).unwrap();
        assert_eq!(sys.len(), 3 * 27);
        let oh = (sys.positions[1] - sys.positions[0]).norm();
        let hh = (sys.positions[1] - sys.positions[2]).norm();
        assert!((oh - WATER_OH).abs() < 1e-9, "O-H {oh}");
        assert!((hh - WATER_HH).abs() < 1e-9, "H-H {hh}");
    }

    #[test]
    fn ion_placement_counts_neutrality_and_separation() {
        let mut sys = place_water_grid(box10(), 0.75).unwrap();
        let n_water = sys.len() / 3;
        let placed = replace_waters_with_ions(&mut sys, 4, 2024, 0.5).unwrap();
        assert_eq!(placed.sodium, 4);
        assert_eq!(placed.chloride, 4);
        assert_eq!(sys.len(), 3 * (n_water - 8) + 8);
        let q: f64 = sys.charges.iter().sum();
        assert!(q.abs() < 1e-12, "net charge {q}");
        let ions: Vec<Vec3> = sys
            .names
            .iter()
            .zip(&sys.positions)
            .filter(|(n, _)| n.as_str() == "NA" || n.as_str() == "CL")
            .map(|(_, p)| *p)
            .collect();
        assert_eq!(ions.len(), 8);
        for (a, pa) in ions.iter().enumerate() {
            for pb in &ions[a + 1..] {
                let d = sys.sim_box.min_image(*pa - *pb).norm();
                assert!(d >= 0.5 - 1e-12, "ions {d} nm apart");
            }
        }
    }

    #[test]
    fn ion_placement_is_deterministic_for_a_seed() {
        let mut a = place_water_grid(box10(), 0.75).unwrap();
        let mut b = place_water_grid(box10(), 0.75).unwrap();
        replace_waters_with_ions(&mut a, 3, 7, 0.2).unwrap();
        replace_waters_with_ions(&mut b, 3, 7, 0.2).unwrap();
        assert_eq!(a.positions, b.positions);
        let mut c = place_water_grid(box10(), 0.75).unwrap();
        replace_waters_with_ions(&mut c, 3, 8, 0.2).unwrap();
        assert_ne!(a.positions, c.positions);
    }

    #[test]
    fn impossible_separation_is_reported() {
        let mut sys = place_water_grid(box10(), 1.0).unwrap();
        assert!(replace_waters_with_ions(&mut sys, 20, 1, 2.0).is_err());
    }
}

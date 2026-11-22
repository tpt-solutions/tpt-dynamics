//! GROMACS `.gro` reader and writer (fixed-width, nm and nm/ps units).
//!
//! Layout: a title line, the atom count, one fixed-width line per atom
//! (`%5d%-5s%5s%5d%8.3f%8.3f%8.3f` plus optional `%8.4f` velocities), and a
//! final box line. The box line holds three values for an orthorhombic box, or nine values
//! (GROMACS order `v1x v2y v3z v1y v1z v2x v2z v3x v3y`) for a triclinic box.
//! The nine-value form must be lower-triangular: `v1y`, `v1z` and `v2z` are zero.
//!
//! Atom and residue numbers are written modulo 100 000 as the format requires.

use std::io::{BufRead, Write};

use tpt_dynamics_core::Vec3;

use crate::error::SimError;
use crate::pbc::SimBox;
use crate::system::{element_mass, guess_element, NewAtom, System};
use crate::textio::{BoundedLines, ReadLimits};

const FORMAT: &str = "GRO";

/// Returns the substring in byte range `start..end` (0-based, exclusive end),
/// or an empty string if the range is out of bounds. Never panics.
fn field(line: &str, start: usize, end: usize) -> &str {
    line.get(start..end.min(line.len())).unwrap_or("")
}

fn parse_f64(text: &str, line_no: usize, what: &str) -> Result<f64, SimError> {
    let v: f64 = text.trim().parse().map_err(|_| {
        SimError::parse(FORMAT, line_no, format!("invalid {what} '{}'", text.trim()))
    })?;
    if !v.is_finite() {
        return Err(SimError::parse(FORMAT, line_no, format!("non-finite {what}")));
    }
    Ok(v)
}

/// Reads a GRO structure from a buffered reader, with the default size limits
/// ([`ReadLimits::default`]).
///
/// Velocities are read when present on every atom line; otherwise they are zero.
///
/// # Errors
/// Returns [`SimError::Parse`] for malformed content (bad counts, numbers,
/// missing or non-orthorhombic box, unknown element, a line or atom count over
/// the limits) and [`SimError::Io`] on read failure.
pub fn read_gro<R: BufRead>(reader: R) -> Result<System, SimError> {
    read_gro_with_limits(reader, ReadLimits::default())
}

/// Reads a GRO structure, rejecting lines longer than `limits.max_line_bytes`
/// and files declaring more than `limits.max_atoms` atoms.
///
/// # Errors
/// As [`read_gro`].
pub fn read_gro_with_limits<R: BufRead>(reader: R, limits: ReadLimits) -> Result<System, SimError> {
    let mut lines = BoundedLines::new(reader, FORMAT, limits.max_line_bytes);
    let mut line_no = 0usize;

    let mut next_line = |line_no: &mut usize| -> Result<String, SimError> {
        *line_no += 1;
        match lines.next() {
            Some(l) => l,
            None => Err(SimError::parse(FORMAT, *line_no, "unexpected end of file")),
        }
    };

    let _title = next_line(&mut line_no)?;
    let count_line = next_line(&mut line_no)?;
    let count: usize = count_line.trim().parse().map_err(|_| {
        SimError::parse(FORMAT, line_no, format!("invalid atom count '{}'", count_line.trim()))
    })?;
    if count == 0 {
        return Err(SimError::parse(FORMAT, line_no, "atom count is zero"));
    }
    if count > limits.max_atoms {
        return Err(SimError::parse(
            FORMAT,
            line_no,
            format!("atom count {count} exceeds the limit of {}", limits.max_atoms),
        ));
    }

    let mut atoms = Vec::with_capacity(count.min(1 << 20));
    for _ in 0..count {
        let line = next_line(&mut line_no)?;
        if line.len() < 44 {
            return Err(SimError::parse(FORMAT, line_no, "atom line shorter than 44 columns"));
        }
        let residue_id: i32 = {
            let t = field(&line, 0, 5).trim();
            if t.is_empty() {
                0
            } else {
                t.parse().map_err(|_| {
                    SimError::parse(FORMAT, line_no, format!("invalid residue number '{t}'"))
                })?
            }
        };
        let residue_name = field(&line, 5, 10).trim().to_string();
        let name = field(&line, 10, 15).trim().to_string();
        let x = parse_f64(field(&line, 20, 28), line_no, "x coordinate")?;
        let y = parse_f64(field(&line, 28, 36), line_no, "y coordinate")?;
        let z = parse_f64(field(&line, 36, 44), line_no, "z coordinate")?;

        let v_fields = [field(&line, 44, 52), field(&line, 52, 60), field(&line, 60, 68)];
        let present = v_fields.iter().filter(|s| !s.trim().is_empty()).count();
        let velocity = match present {
            0 => Vec3::ZERO,
            3 => Vec3::new(
                parse_f64(v_fields[0], line_no, "vx")?,
                parse_f64(v_fields[1], line_no, "vy")?,
                parse_f64(v_fields[2], line_no, "vz")?,
            ),
            _ => return Err(SimError::parse(FORMAT, line_no, "partial velocity fields")),
        };

        let element = guess_element(&name, &residue_name).ok_or_else(|| {
            SimError::parse(FORMAT, line_no, format!("cannot guess element for '{name}'"))
        })?;
        let mass = element_mass(&element).ok_or_else(|| {
            SimError::parse(
                FORMAT,
                line_no,
                format!("unknown element '{element}' for '{name}', cannot assign mass"),
            )
        })?;
        atoms.push(NewAtom {
            name,
            residue_name,
            residue_id,
            element,
            position: Vec3::new(x, y, z),
            velocity,
            mass,
            charge: 0.0,
            atom_type: 0,
        });
    }

    let box_line = next_line(&mut line_no)?;
    let values: Vec<f64> = box_line
        .split_whitespace()
        .map(|t| parse_f64(t, line_no, "box value"))
        .collect::<Result<_, _>>()?;
    let sim_box = parse_box(&values, line_no)?;

    let mut system = System::new(sim_box);
    for atom in atoms {
        system.push_atom(atom);
    }
    Ok(system)
}

/// Parses the GRO box line: three values (orthorhombic) or nine values in
/// GROMACS order `v1x v2y v3z v1y v1z v2x v2z v3x v3y`.
///
/// The nine-value form is lower-triangular: v1 = (lx, 0, 0), v2 = (xy, ly, 0),
/// v3 = (xz, yz, lz). So `v1y`, `v1z` and `v2z` must be zero.
fn parse_box(values: &[f64], line_no: usize) -> Result<SimBox, SimError> {
    let sim_box = match values {
        [lx, ly, lz] => SimBox::new(*lx, *ly, *lz),
        [v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y] => {
            if *v1y != 0.0 || *v1z != 0.0 || *v2z != 0.0 {
                return Err(SimError::parse(
                    FORMAT,
                    line_no,
                    "box vectors must be lower-triangular (v1 along x, v2 in the xy plane)",
                ));
            }
            SimBox::triclinic(*v1x, *v2y, *v3z, *v2x, *v3x, *v3y)
        }
        _ => return Err(SimError::parse(FORMAT, line_no, "box line must have 3 or 9 values")),
    };
    sim_box.map_err(|e| SimError::parse(FORMAT, line_no, e.to_string()))
}
fn fixed(value: f64, width: usize, decimals: usize) -> Result<String, SimError> {
    let s = format!("{value:>width$.decimals$}");
    if s.len() != width {
        return Err(SimError::Invalid(format!("value {value} does not fit in {width} columns")));
    }
    Ok(s)
}

/// Writes a system as a GRO file with positions (nm, 3 decimals), velocities
/// (nm/ps, 4 decimals) and an orthorhombic box line.
///
/// # Errors
/// Returns [`SimError::Invalid`] if a name is too long for its field or a value
/// does not fit its column width, and [`SimError::Io`] on write failure.
pub fn write_gro<W: Write>(mut writer: W, system: &System, title: &str) -> Result<(), SimError> {
    writeln!(writer, "{}", title.replace(['\n', '\r'], " "))?;
    writeln!(writer, "{:>5}", system.len())?;
    for i in 0..system.len() {
        let name = &system.names[i];
        let res = &system.residue_names[i];
        if name.len() > 5 || res.len() > 5 {
            return Err(SimError::Invalid(format!("atom {i}: name longer than 5 characters")));
        }
        let p = system.positions[i];
        let v = system.velocities[i];
        let resnum = system.residue_ids[i].rem_euclid(100_000);
        let serial = (i + 1) % 100_000;
        writeln!(
            writer,
            "{resnum:>5}{res:<5}{name:>5}{serial:>5}{}{}{}{}{}{}",
            fixed(p.x, 8, 3)?,
            fixed(p.y, 8, 3)?,
            fixed(p.z, 8, 3)?,
            fixed(v.x, 8, 4)?,
            fixed(v.y, 8, 4)?,
            fixed(v.z, 8, 4)?,
        )?;
    }
    let b = &system.sim_box;
    if b.is_orthorhombic() {
        let l = b.lengths;
        writeln!(writer, "{:>10.5}{:>10.5}{:>10.5}", l.x, l.y, l.z)?;
    } else {
        // Nine values in GROMACS order: v1x v2y v3z v1y v1z v2x v2z v3x v3y.
        let l = b.lengths;
        let [xy, xz, yz] = b.tilt;
        let values = [l.x, l.y, l.z, 0.0, 0.0, xy, 0.0, xz, yz];
        for v in values {
            write!(writer, "{v:>10.5}")?;
        }
        writeln!(writer)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::NewAtom;

    fn sample() -> System {
        let mut s = System::new(SimBox::new(2.5, 3.0, 3.5).unwrap());
        let specs = [("OW", "SOL", "O", 1), ("HW1", "SOL", "H", 1), ("NA", "NA", "Na", 2)];
        for (i, (name, res, el, rid)) in specs.iter().enumerate() {
            s.push_atom(NewAtom {
                name: (*name).into(),
                residue_name: (*res).into(),
                residue_id: *rid,
                element: (*el).into(),
                position: Vec3::new(0.1 + i as f64 * 0.7, 1.234 + i as f64, 2.0 - 0.25 * i as f64),
                velocity: Vec3::new(0.1234 * (i as f64 + 1.0), -0.5, 0.0009),
                mass: element_mass(el).unwrap(),
                charge: 0.0,
                atom_type: 0,
            });
        }
        s
    }

    #[test]
    fn round_trip_preserves_state() {
        let s = sample();
        let mut buf = Vec::new();
        write_gro(&mut buf, &s, "test system").unwrap();
        let back = read_gro(&buf[..]).unwrap();
        assert_eq!(back.len(), 3);
        assert_eq!(back.names, s.names);
        assert_eq!(back.residue_names, s.residue_names);
        assert_eq!(back.residue_ids, s.residue_ids);
        assert_eq!(back.elements, s.elements);
        for i in 0..s.len() {
            assert!((back.positions[i] - s.positions[i]).norm() < 5e-4);
            assert!((back.velocities[i] - s.velocities[i]).norm() < 5e-5);
        }
        assert!((back.sim_box.lengths.y - 3.0).abs() < 1e-9);
    }

    #[test]
    fn reads_without_velocities_and_nine_value_box() {
        let text = "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   2.00000   2.00000   2.00000   0.00000   0.00000   0.00000   0.00000   0.00000   0.00000\n";
        let s = read_gro(text.as_bytes()).unwrap();
        assert_eq!(s.len(), 1);
        assert_eq!(s.velocities[0], Vec3::ZERO);
        assert!((s.positions[0].x - 0.126).abs() < 1e-12);
        assert_eq!(s.elements[0], "O");
    }

    #[test]
    fn malformed_input_is_an_error_not_a_panic() {
        let cases = [
            "",
            "title\n",
            "title\nabc\n",
            "title\n    0\n   1.0 1.0 1.0\n",
            "title\n    1\n    1SOL     OW    1   0.126   1.624\n   2.0 2.0 2.0\n",
            "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n",
            "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   2.0 2.0\n",
            "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   2.0 1.0 2.0 0.5 0 0 0 0 0\n",
            "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679   0.1\n   2.0 2.0 2.0\n",
            "title\n    1\n    1SOL     XX    1   0.126   1.624   1.679\n   2.0 2.0 2.0\n",
            "title\n    1\n    1SOL     OW    1   0.126   nan   1.679\n   2.0 2.0 2.0\n",
            "title\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   0.0 2.0 2.0\n",
        ];
        for case in cases {
            assert!(read_gro(case.as_bytes()).is_err(), "expected error for {case:?}");
        }
    }

    #[test]
    fn every_truncation_of_valid_input_parses_without_panic() {
        let mut buf = Vec::new();
        write_gro(&mut buf, &sample(), "t").unwrap();
        for n in 0..buf.len() {
            let _ = read_gro(&buf[..n]);
        }
    }
}

#[cfg(test)]
mod limit_tests {
    use super::*;
    use crate::textio::ReadLimits;

    #[test]
    fn overlong_line_is_rejected() {
        let mut text = b"title\n    1\n".to_vec();
        text.extend(std::iter::repeat_n(b' ', 2 << 20));
        text.extend_from_slice(b"\n   2.0 2.0 2.0\n");
        let err = read_gro(&text[..]).unwrap_err().to_string();
        assert!(err.contains("line 3") && err.contains("longer than"), "{err}");
    }

    #[test]
    fn declared_atom_count_over_limit_is_rejected_before_reading_atoms() {
        let text = "title\n   10\n";
        let limits = ReadLimits { max_line_bytes: 1 << 20, max_atoms: 9 };
        let err = read_gro_with_limits(text.as_bytes(), limits).unwrap_err().to_string();
        assert!(err.contains("exceeds the limit of 9"), "{err}");
        let huge = "title\n99999999999\n";
        assert!(read_gro(huge.as_bytes()).is_err());
    }
}

#[cfg(test)]
mod triclinic_tests {
    use super::*;
    use crate::system::NewAtom;

    fn tri_system(tilt: [f64; 3]) -> System {
        let [xy, xz, yz] = tilt;
        let mut s = System::new(SimBox::triclinic(3.0, 2.5, 4.0, xy, xz, yz).unwrap());
        s.push_atom(NewAtom {
            name: "OW".into(),
            residue_name: "SOL".into(),
            residue_id: 1,
            element: "O".into(),
            position: Vec3::new(0.5, 1.0, 2.0),
            velocity: Vec3::ZERO,
            mass: element_mass("O").unwrap(),
            charge: 0.0,
            atom_type: 0,
        });
        s
    }

    #[test]
    fn triclinic_box_round_trips_through_nine_values() {
        let s = tri_system([1.0, -0.7, 0.4]);
        let mut buf = Vec::new();
        write_gro(&mut buf, &s, "tri").unwrap();
        let text = String::from_utf8(buf.clone()).unwrap();
        assert_eq!(text.lines().last().unwrap().split_whitespace().count(), 9);
        let back = read_gro(&buf[..]).unwrap();
        assert!((back.sim_box.lengths - s.sim_box.lengths).norm() < 1e-5);
        for (a, b) in back.sim_box.tilt.iter().zip(&s.sim_box.tilt) {
            assert!((a - b).abs() < 1e-5, "tilt {a} vs {b}");
        }
        assert!(!back.sim_box.is_orthorhombic());
    }

    #[test]
    fn orthorhombic_box_is_still_written_with_three_values() {
        let mut buf = Vec::new();
        write_gro(&mut buf, &sample_ortho(), "o").unwrap();
        let text = String::from_utf8(buf).unwrap();
        assert_eq!(text.lines().last().unwrap().split_whitespace().count(), 3);
    }

    fn sample_ortho() -> System {
        System::new(SimBox::new(2.0, 2.0, 2.0).unwrap())
    }

    #[test]
    fn rejects_nine_value_box_that_is_not_lower_triangular() {
        // v1y (fourth value) nonzero: first vector not along x.
        let text = "t\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   2.0 2.0 2.0 0.1 0.0 0.0 0.0 0.0 0.0\n";
        let err = read_gro(text.as_bytes()).unwrap_err().to_string();
        assert!(err.contains("lower-triangular"), "{err}");
        // v2z (seventh value) nonzero: second vector out of the xy plane.
        let text = "t\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   2.0 2.0 2.0 0.0 0.0 0.0 0.1 0.0 0.0\n";
        assert!(read_gro(text.as_bytes()).is_err());
    }

    #[test]
    fn reads_triclinic_box_values_in_gromacs_order() {
        // v2x = 0.3 (xy), v3x = -0.2 (xz), v3y = 0.1 (yz).
        let text = "t\n    1\n    1SOL     OW    1   0.126   1.624   1.679\n   2.0 2.5 4.0 0.0 0.0 0.3 0.0 -0.2 0.1\n";
        let s = read_gro(text.as_bytes()).unwrap();
        assert_eq!(s.sim_box.tilt, [0.3, -0.2, 0.1]);
        assert_eq!(s.sim_box.lengths, Vec3::new(2.0, 2.5, 4.0));
    }
}

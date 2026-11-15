//! PDB reader and writer (fixed-column format, coordinates in Å).
//!
//! Reading:
//! - `CRYST1` gives the orthorhombic box (required; non-90° angles are rejected).
//! - `ATOM` and `HETATM` records are both read as atoms. Only the first model is
//!   read (reading stops at `ENDMDL` or `END`).
//! - Element comes from columns 77–78 or is guessed from the atom name.
//! - Masses are looked up from the element; unknown elements are an error.
//!
//! Writing emits a `CRYST1` record followed by `ATOM` records (no `HETATM`
//! distinction is kept) and `END`. Serial numbers wrap at 100 000 and residue
//! numbers at 10 000, as the fixed-width format requires.

use std::io::{BufRead, Write};

use tpt_dynamics_core::{units::angstrom_to_nm, Vec3};

use crate::error::SimError;
use crate::pbc::SimBox;
use crate::system::{element_mass, guess_element, normalise_element, NewAtom, System};
use crate::textio::{BoundedLines, ReadLimits};

const FORMAT: &str = "PDB";

/// Returns the trimmed text in 1-based inclusive columns `start..=end`,
/// or an empty string if the line is too short. Never panics.
fn col(line: &str, start: usize, end: usize) -> &str {
    line.get(start - 1..end.min(line.len())).unwrap_or("").trim()
}

fn parse_f64(text: &str, line_no: usize, what: &str) -> Result<f64, SimError> {
    let v: f64 = text
        .parse()
        .map_err(|_| SimError::parse(FORMAT, line_no, format!("invalid {what} '{text}'")))?;
    if !v.is_finite() {
        return Err(SimError::parse(FORMAT, line_no, format!("non-finite {what}")));
    }
    Ok(v)
}

/// Reads a PDB structure from a buffered reader, with the default size limits
/// ([`ReadLimits::default`]).
///
/// # Errors
/// Returns [`SimError::Parse`] for malformed records (missing box, bad numbers,
/// unknown element, no atoms, a line or atom count over the limits) and
/// [`SimError::Io`] for read failures.
pub fn read_pdb<R: BufRead>(reader: R) -> Result<System, SimError> {
    read_pdb_with_limits(reader, ReadLimits::default())
}

/// Reads a PDB structure, rejecting lines longer than `limits.max_line_bytes`
/// and files with more than `limits.max_atoms` atoms.
///
/// # Errors
/// As [`read_pdb`].
pub fn read_pdb_with_limits<R: BufRead>(reader: R, limits: ReadLimits) -> Result<System, SimError> {
    let mut sim_box: Option<SimBox> = None;
    let mut atoms: Vec<NewAtom> = Vec::new();

    for (idx, line) in BoundedLines::new(reader, FORMAT, limits.max_line_bytes).enumerate() {
        let line_no = idx + 1;
        let line = line?;
        let record = col(&line, 1, 6);
        match record {
            "CRYST1" => {
                let a = parse_f64(col(&line, 7, 15), line_no, "box a")?;
                let b = parse_f64(col(&line, 16, 24), line_no, "box b")?;
                let c = parse_f64(col(&line, 25, 33), line_no, "box c")?;
                let alpha = parse_f64(col(&line, 34, 40), line_no, "alpha")?;
                let beta = parse_f64(col(&line, 41, 47), line_no, "beta")?;
                let gamma = parse_f64(col(&line, 48, 54), line_no, "gamma")?;
                let right_angles = [alpha, beta, gamma].iter().all(|x| (x - 90.0).abs() < 1e-3);
                if !right_angles {
                    return Err(SimError::parse(
                        FORMAT,
                        line_no,
                        "non-orthorhombic CRYST1 cell is not supported",
                    ));
                }
                let bx = SimBox::new(angstrom_to_nm(a), angstrom_to_nm(b), angstrom_to_nm(c))
                    .map_err(|e| SimError::parse(FORMAT, line_no, e.to_string()))?;
                sim_box = Some(bx);
            }
            "ATOM" | "HETATM" => {
                if line.len() < 54 {
                    return Err(SimError::parse(
                        FORMAT,
                        line_no,
                        "atom record shorter than 54 columns",
                    ));
                }
                let x = parse_f64(col(&line, 31, 38), line_no, "x coordinate")?;
                let y = parse_f64(col(&line, 39, 46), line_no, "y coordinate")?;
                let z = parse_f64(col(&line, 47, 54), line_no, "z coordinate")?;
                let name = col(&line, 13, 16).to_string();
                let residue_name = col(&line, 18, 20).to_string();
                let res_text = col(&line, 23, 26);
                let residue_id = if res_text.is_empty() {
                    0
                } else {
                    res_text.parse::<i32>().map_err(|_| {
                        SimError::parse(
                            FORMAT,
                            line_no,
                            format!("invalid residue number '{res_text}'"),
                        )
                    })?
                };
                let elem_text = col(&line, 77, 78);
                let element = if elem_text.is_empty() {
                    guess_element(&name, &residue_name).unwrap_or_default()
                } else {
                    normalise_element(elem_text)
                };
                let mass = element_mass(&element).ok_or_else(|| {
                    SimError::parse(
                        FORMAT,
                        line_no,
                        format!("unknown element '{element}', cannot assign mass"),
                    )
                })?;
                if atoms.len() >= limits.max_atoms {
                    return Err(SimError::parse(
                        FORMAT,
                        line_no,
                        format!("file has more than {} atoms", limits.max_atoms),
                    ));
                }
                atoms.push(NewAtom {
                    name,
                    residue_name,
                    residue_id,
                    element,
                    position: Vec3::new(angstrom_to_nm(x), angstrom_to_nm(y), angstrom_to_nm(z)),
                    velocity: Vec3::ZERO,
                    mass,
                    charge: 0.0,
                    atom_type: 0,
                });
            }
            "ENDMDL" | "END" if !atoms.is_empty() => break,
            _ => {}
        }
    }

    let sim_box = sim_box
        .ok_or_else(|| SimError::parse(FORMAT, 0, "missing CRYST1 record (box required)"))?;
    if atoms.is_empty() {
        return Err(SimError::parse(FORMAT, 0, "no ATOM or HETATM records found"));
    }
    let mut system = System::new(sim_box);
    for atom in atoms {
        system.push_atom(atom);
    }
    Ok(system)
}

fn checked(field: &str, value: f64, width: usize, decimals: usize) -> Result<String, SimError> {
    let s = format!("{value:>width$.decimals$}");
    if s.len() != width {
        return Err(SimError::Invalid(format!("{field} {value} does not fit in {width} columns")));
    }
    Ok(s)
}

/// Writes a system as a PDB file (one model, `CRYST1` + `ATOM` records).
///
/// Positions are written in Å with three decimals. Velocities are not stored
/// in PDB and are dropped.
///
/// # Errors
/// Returns [`SimError::Invalid`] if a name or coordinate does not fit its
/// fixed-width field, and [`SimError::Io`] on write failure.
pub fn write_pdb<W: Write>(mut writer: W, system: &System) -> Result<(), SimError> {
    let l = system.sim_box.lengths;
    writeln!(
        writer,
        "CRYST1{:>9.3}{:>9.3}{:>9.3}{:>7.2}{:>7.2}{:>7.2} {:<11}{:>4}",
        l.x * 10.0,
        l.y * 10.0,
        l.z * 10.0,
        90.0,
        90.0,
        90.0,
        "P 1",
        1
    )?;
    for i in 0..system.len() {
        let name = &system.names[i];
        let res = &system.residue_names[i];
        let element = &system.elements[i];
        if name.len() > 4 || res.len() > 3 || element.len() > 2 {
            return Err(SimError::Invalid(format!(
                "atom {i}: name, residue name or element too long for PDB fields"
            )));
        }
        // PDB convention: one-letter element symbols with short names start in column 14.
        let name_field = if element.len() == 1 && name.len() < 4 {
            format!(" {name:<3}")
        } else {
            format!("{name:<4}")
        };
        let p = system.positions[i];
        let x = checked("x", p.x * 10.0, 8, 3)?;
        let y = checked("y", p.y * 10.0, 8, 3)?;
        let z = checked("z", p.z * 10.0, 8, 3)?;
        let serial = (i + 1) % 100_000;
        let resnum = system.residue_ids[i].rem_euclid(10_000);
        writeln!(
            writer,
            "ATOM  {serial:>5} {name_field} {res:>3} A{resnum:>4}    {x}{y}{z}{:>6.2}{:>6.2}          {element:>2}",
            1.0, 0.0
        )?;
    }
    writeln!(writer, "END")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::system::NewAtom;

    fn sample() -> System {
        let mut s = System::new(SimBox::new(3.0, 4.0, 5.0).unwrap());
        for (i, (name, res, el)) in
            [("N", "ALA", "N"), ("CA", "ALA", "C"), ("OW", "SOL", "O"), ("CL", "CL", "Cl")]
                .iter()
                .enumerate()
        {
            s.push_atom(NewAtom {
                name: (*name).into(),
                residue_name: (*res).into(),
                residue_id: i32::try_from(i + 1).unwrap(),
                element: (*el).into(),
                position: Vec3::new(0.1 * i as f64 + 0.123, 1.0 + 0.5 * i as f64, 2.25),
                velocity: Vec3::ZERO,
                mass: element_mass(el).unwrap(),
                charge: 0.0,
                atom_type: 0,
            });
        }
        s
    }

    #[test]
    fn round_trip_preserves_structure() {
        let s = sample();
        let mut buf = Vec::new();
        write_pdb(&mut buf, &s).unwrap();
        let text = String::from_utf8(buf.clone()).unwrap();
        assert!(text.lines().next().unwrap().starts_with("CRYST1"));
        let back = read_pdb(&buf[..]).unwrap();
        assert_eq!(back.len(), s.len());
        assert_eq!(back.names, s.names);
        assert_eq!(back.residue_names, s.residue_names);
        assert_eq!(back.residue_ids, s.residue_ids);
        assert_eq!(back.elements, s.elements);
        for (a, b) in back.positions.iter().zip(&s.positions) {
            assert!((*a - *b).norm() < 5e-4, "{a:?} vs {b:?}");
        }
        assert!((back.sim_box.lengths.z - 5.0).abs() < 1e-9);
    }

    #[test]
    fn reads_fixed_column_record() {
        let text = "CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n\
ATOM      1  OW  SOL A   1      10.000  20.000  30.000  1.00  0.00           O\n\
HETATM    2 CL   CL B   2       1.500   2.500   3.500  1.00  0.00           CL\n\
END\n";
        let s = read_pdb(text.as_bytes()).unwrap();
        assert_eq!(s.len(), 2);
        assert_eq!(s.names[1], "CL");
        assert_eq!(s.elements[1], "Cl");
        assert!((s.positions[0].x - 1.0).abs() < 1e-12);
        assert!((s.positions[1].z - 0.35).abs() < 1e-12);
    }

    #[test]
    fn malformed_input_is_an_error_not_a_panic() {
        let cases = [
            "",
            "ATOM      1  N   ALA A   1      1.000   2.000   3.000\n",
            "CRYST1   30.000   40.000   50.000  80.00  90.00  90.00 P 1           1\n",
            "CRYST1   30.000   40.000\n",
            "CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n\
ATOM      1  N   ALA A   1      abc     2.000   3.000  1.00  0.00           N\n",
            "CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n\
ATOM      1  XX  ALA A   1      1.000   2.000   3.000  1.00  0.00          XX\n",
            "CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n\
ATOM      1  N   ALA A  zz      1.000   2.000   3.000  1.00  0.00           N\n",
            "CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n",
            "CRYST1   0.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n\
ATOM      1  N   ALA A   1      1.000   2.000   3.000  1.00  0.00           N\n",
        ];
        for case in cases {
            assert!(read_pdb(case.as_bytes()).is_err(), "expected error for {case:?}");
        }
    }

    #[test]
    fn every_truncation_of_valid_input_parses_without_panic() {
        let mut buf = Vec::new();
        write_pdb(&mut buf, &sample()).unwrap();
        for n in 0..buf.len() {
            let _ = read_pdb(&buf[..n]);
        }
    }

    #[test]
    fn overflowing_coordinate_is_rejected() {
        let mut s = sample();
        s.positions[0] = Vec3::new(2000.0, 0.0, 0.0);
        let mut buf = Vec::new();
        assert!(write_pdb(&mut buf, &s).is_err());
    }
}

#[cfg(test)]
mod limit_tests {
    use super::*;
    use crate::textio::ReadLimits;

    const HEADER: &str = "CRYST1   30.000   40.000   50.000  90.00  90.00  90.00 P 1           1\n";
    const ATOM: &str =
        "ATOM      1  OW  SOL A   1      10.000  20.000  30.000  1.00  0.00           O\n";

    #[test]
    fn overlong_line_is_rejected_with_its_line_number() {
        // A 2 MiB comment line with no newline inside it (security finding S-1).
        let mut text = HEADER.as_bytes().to_vec();
        text.extend_from_slice(b"REMARK ");
        text.extend(std::iter::repeat_n(b'x', 2 << 20));
        text.extend_from_slice(ATOM.as_bytes());
        let err = read_pdb(&text[..]).unwrap_err().to_string();
        assert!(err.contains("line 2") && err.contains("longer than"), "{err}");
    }

    #[test]
    fn atom_count_limit_is_enforced() {
        let text = format!("{HEADER}{ATOM}{ATOM}{ATOM}END\n");
        let limits = ReadLimits { max_line_bytes: 1 << 20, max_atoms: 2 };
        let err = read_pdb_with_limits(text.as_bytes(), limits).unwrap_err().to_string();
        assert!(err.contains("more than 2 atoms"), "{err}");
        let ok = ReadLimits { max_line_bytes: 1 << 20, max_atoms: 3 };
        assert_eq!(read_pdb_with_limits(text.as_bytes(), ok).unwrap().len(), 3);
    }
}

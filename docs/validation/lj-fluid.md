# Lennard-Jones fluid: equation-of-state check

Status: **simulation run; literature comparison inconclusive**. The NIST table is
now read and checked, but its T* = 1.00 entries near ρ* = 0.8 are not confirmed
as a pass/fail reference, and the raw data files could not be reconciled with the
table (see "Reference status"). Do not treat this file as a validation pass.

## Setup

- Code: `tpt-dynamics-sim/examples/lj_eos_check.rs`, release build.
- Potential: full 12-6 Lennard-Jones, truncated (unshifted) at 3σ, with
  analytic long-range energy and pressure corrections added afterwards. This is
  the convention "3σ + standard long range corrections" used by the NIST TMMC
  tables.
- State point: reduced density ρ* = 0.8, reduced temperature T* = 1.0.
- System: 864 atoms, FCC start, cubic periodic box.
- Dynamics: velocity Verlet, dt = 0.002 (reduced), stochastic velocity
  rescaling to T* = 1.0 (τ = 0.5). Equilibration 5000 steps, production
  10 blocks of 1000 steps, sampled every 10 steps.

## Result

| Quantity | Simulated (± block SEM, 10 blocks) |
|---|---|
| P* | 1.0685 ± 0.0147 |
| U*/N (with tail) | −5.5216 ± 0.0037 |

Tail corrections applied: U_tail/N = −0.2481, P_tail = −0.3968 (analytic; the
formulas were rechecked in 2026-11-20 against the standard expressions, and the
values agree).

Run: 44 s wall time for the whole example, release build, single run.

## Reference status

Source: NIST LJ benchmark, EOS-TMMC page, <https://mmlapps.nist.gov/srs/LJ_PURE/eostmmc.htm>,
"Truncation: 3σ + standard long range corrections", V/σ³ = 512. Read on 2026-11-20 with
the web-fetch agent, which returns summaries. Exact numbers below are verbatim from the page
text, but the reading was done through a summarising tool, so re-check any row against the
raw HTML before using it as a pass/fail value.

- The page lists per-temperature blocks headed `kBT/ε= X`, with columns
  `ρσ3  pσ3/ε  +/-`. The ± column is one standard deviation from five independent
  simulations (67 % confidence), per the page remarks.
- The T* = 1.00 block has 410 rows, at densities that are multiples of 1/512. The
  rows bracketing ρ* = 0.8 are:
  - ρσ³ = 0.798828, pσ³/ε = 0.988636 (±0.00655 per the earlier reading; re-check)
  - ρσ³ = 0.800781, pσ³/ε = 1.01472 (±0.00396)
  Linear interpolation to ρσ³ = 0.8 gives **P* ≈ 1.004**. This is an interpolation,
  not a printed value.
- The first row of the block, ρσ³ = 1.95313E-03 with pσ³/ε = 6.56814E-04, is **not**
  consistent with the ideal-gas limit (P* ≈ ρ* ≈ 1.95E-03 at this density). Either the
  printed pressure column is wrong at low density or the columns are mislabelled.
  This casts doubt on the column assignment for the whole block.
- The earlier candidate in this file (ρσ³ = 0.80078, pσ³/ε = 0.98864) paired the
  wrong density with the pressure. The 0.98864 value belongs to ρσ³ = 0.798828.
  The earlier "about 8 % higher" statement is withdrawn.
- No U*/N table is printed on the page. The page links raw TMMC data per temperature
  (for T* = 1.00: `lj.t100.tar.gz`, containing five runs of ln Π(N) and average total
  potential energy per N).

### Raw TMMC data (lj.t100.tar.gz): not yet usable

The raw files were downloaded and extracted into a scratch directory and analysed.
Results:

- The ideal-gas check at N = 1 matches to 5×10⁻⁴ when ln Π(N) is read as
  ln Q_N + N ln z_sim + const, with ln z_sim = −3.823078 (the value on the page).
- At larger N, ln Q_N − (ideal gas) grows as about 5 N²/V. That is far larger than
  the expected long-range interaction term (about 0.31 N²/V for these tail
  corrections), so the raw files carry a convention I have not identified.
- Consequently, pressures obtained from the raw files do not pass a dilute-gas check
  at ρ* = 0.1 (P* = 0.025 against 0.1 expected). At ρ* ≈ 0.8 the raw-file values
  (P* ≈ 1.017, U*/N ≈ −5.53 before tail correction) are not trusted either.

Until the TMMC convention for these files is confirmed from the source paper (Johnson,
Zollweg & Gubbins 1993, not yet read), the raw data is not a reference.

## Conclusion

The simulation is stable and gives a consistent, physically plausible state
(P* ≈ 1.07, U*/N ≈ −5.52 with tail corrections). The comparison is **not
established either way**:

- Against the interpolated table value (P* ≈ 1.004), the simulated P* is about 6 %
  higher. That gap is larger than the stated block SEM, but the block SEM comes from
  only 10 blocks and may be underestimated, and the reference values carry their own
  uncertainty. Treat this as an open discrepancy to investigate, not a confirmed error.
- The table's low-density row fails the ideal-gas limit, so the column assignment in
  that block is unresolved.

Next steps to close this: (1) read Johnson et al. 1993 for the TMMC convention and the
exact T* = 1.00 table, (2) rerun the simulation with longer production (at least 50
blocks) and check block-SEM convergence, (3) confirm the energy convention (truncated
vs. with tail) before comparing U*/N.

## Reproduce

```text
cargo run --release -p tpt-dynamics-sim --example lj_eos_check
```

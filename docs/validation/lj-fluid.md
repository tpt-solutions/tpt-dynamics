# Lennard-Jones fluid: equation-of-state check

Status: **simulation run; literature comparison not completed** (no verified
reference value; see "Reference status").

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

Tail corrections applied: U_tail/N = −0.2481, P_tail = −0.3968 (analytic,
checked by hand against the standard formulas).

Run: 44 s wall time for the whole example, release build, single run.

## Reference status

- Johnson, Zollweg & Gubbins, Mol. Phys. 78, 591 (1993), is the target EOS.
  Its tables were not readable through the tools available here (the PDF
  stream could not be decoded), so **no P* or U* value was taken from it**.
- NIST LJ benchmark, EOS-TMMC page (`mmlapps.nist.gov/srs/LJ_PURE/eostmmc.htm`,
  "3σ + standard long range corrections", N = 512). The extracted text has an
  unlabelled block that could be the T* = 1.00 table, with a row at
  ρσ³ = 0.80078 and pσ³/ε = 0.98864 (±0.0065). The block's column order and
  its T* label could not be confirmed, so **this row is a candidate only**.
  Against that candidate, the simulated P* is about 8 % higher. The deviation is
  not established as a real discrepancy, because the reference is uncertain.

## Conclusion

The simulation is stable and gives a consistent, physically plausible state
(P* ≈ 1.07, U*/N ≈ −5.52 with tail corrections). The comparison against a
published equation of state is **inconclusive** until a verified table value at
T* = 1.0 and ρ* = 0.8 (or a nearby state point, with the same truncation) is
read from the source. Do not treat this file as a validation pass.

## Reproduce

```text
cargo run --release -p tpt-dynamics-sim --example lj_eos_check
```

Water box (512 rigid three-site molecules)
==========================================

This template is a scaffold. It does NOT include a force field, because tpt-dynamics
does not bundle any third-party parameter set (see docs/forcefield.md). Supply your
own file, water.toml, in the project's ffparams format. Whoever loads a parameter set
must follow that set's licence and citation terms.

1. Create water.toml next to simulate.toml. It must define:
   - atom types for oxygen and hydrogen (mass, sigma and epsilon in nm and kJ/mol),
     with bond and angle types for the O-H bond and the H-O-H angle,
   - one molecule template named SOL with three atoms in this order: O, H1, H2.
     Give each atom a charge and reference coordinates. Bonds are [[0, 1], [0, 2]].

   The format is described in the module documentation of
   tpt-dynamics-sim/src/ffparams.rs. Check that the SOL geometry matches the water
   model named in [constraints] (water = "tip3p" by default), because SETTLE holds
   the O-H distance and H-O-H angle of that model.

2. Check the configuration without running it:

       tptdyn validate simulate.toml

   The summary should report 1536 atoms (512 molecules x 3).

3. Run it, with an HTML report:

       tptdyn simulate simulate.toml --report

Without water.toml, validate and simulate stop with an input error that names the
missing file. The box is 2.48 nm on each side. The cutoff (0.9 nm) and skin (0.1 nm)
must stay below half the box, which the file already satisfies.

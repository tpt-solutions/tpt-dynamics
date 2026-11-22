Lennard-Jones fluid project (created by `tptdyn init lj-fluid`)

  tptdyn validate simulate.toml          check the configuration, no run
  tptdyn simulate simulate.toml --report run the simulation and write report.html
  tptdyn analyze --trajectory traj.bin --out rmsd.csv
                                         per-frame RMSD and radius of gyration

Outputs: energy.csv (energy log), traj.bin (trajectory), report.html (report).
Edit simulate.toml to change the lattice, temperature, thermostat or run length.

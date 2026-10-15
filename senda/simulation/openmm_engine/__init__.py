"""
senda.simulation.openmm_engine
OpenMM engine for the classical MM stages (01_min, 02_heat, 03_equil, 04_NVT).

Selected with ``amber_simulator.engine: openmm``. The system is still built by
tleap (00_prep is shared with the AMBER engine) and loaded from the same
structure(_HMR).parm7, so force-field parameters and atom order are identical
between engines. Outputs keep the AMBER layout (.rst7 restarts, AMBER NetCDF
.nc trajectories) so senda-analyse and the QM/MM stages consume them unchanged.

Modules
-------
  restraints -- positional and NMR (DISANG &rst) restraints as OpenMM forces
  system     -- System/Context construction from parm7
  io         -- rst7/XML state I/O, AMBER NetCDF writer, reporters
  stages     -- min / heat / equil / prod drivers
  cli        -- senda-omm entry point (run on the compute node)

The restraint builders are engine-agnostic OpenMM forces, so they can be reused
by the ASH QM/MM stages (ASH's OpenMMTheory exposes the underlying System).
"""

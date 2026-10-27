"""
senda.simulation.openmm_engine.system
Build OpenMM Systems/Simulations from the tleap parm7 with AMBER-equivalent settings.

AMBER -> OpenMM mapping
  PME, cut=<A>, vdW tail correction  -> PME, nonbondedCutoff, dispersion correction
                                        (no switching function, as in AMBER)
  ntc=2/ntf=2 (SHAKE H-bonds, SETTLE) -> constraints=HBonds, rigidWater=True
  ntc=1 (minimisation)                -> constraints=None, rigidWater=False
  ntt=3, gamma_ln                     -> LangevinMiddleIntegrator(T, gamma, dt)
  ntp=1, barostat=1|2                 -> MonteCarloBarostat (no Berendsen in OpenMM)
  HMR (hmassrepartition)              -> masses read from structure_HMR.parm7
"""

from __future__ import annotations

from pathlib import Path

import openmm as mm
from openmm import app, unit

from .restraints import nmr_forces, parse_disang, positional_force, read_mask_atoms


def make_system(params: dict, replica_dir: Path, *, cutoff_A: float,
                constrained: bool = True) -> tuple[app.AmberPrmtopFile, mm.System]:
    prmtop = app.AmberPrmtopFile(str(replica_dir / "00_prep" / params["topology"]))
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=cutoff_A * 0.1 * unit.nanometer,
        constraints=app.HBonds if constrained else None,
        rigidWater=constrained,
        ewaldErrorTolerance=float(params["ewald_error_tolerance"]),
    )
    for f in system.getForces():
        if isinstance(f, mm.NonbondedForce):
            f.setUseDispersionCorrection(True)
            f.setUseSwitchingFunction(False)
    if params.get("nmr"):
        for force in nmr_forces(parse_disang(replica_dir / params["nmr"])):
            system.addForce(force)
    return prmtop, system


def add_positional(system: mm.System, params: dict, replica_dir: Path,
                   ref_positions_nm, weight_kcal: float) -> mm.CustomExternalForce:
    indices = read_mask_atoms(replica_dir / params["positional_atoms"])
    force = positional_force(indices, ref_positions_nm, weight_kcal)
    system.addForce(force)
    return force


def add_barostat(system: mm.System, params: dict, temperature_K: float) -> mm.MonteCarloBarostat:
    baro = mm.MonteCarloBarostat(float(params["pressure_bar"]) * unit.bar,
                                 temperature_K * unit.kelvin,
                                 int(params["barostat_interval"]))
    baro.setRandomNumberSeed(int(params["seed"]) % 2**31)
    system.addForce(baro)
    return baro


def platform_and_properties(params: dict) -> tuple[mm.Platform, dict]:
    platform = mm.Platform.getPlatformByName(params["platform"])
    props = {}
    if params["platform"] in ("CUDA", "OpenCL", "HIP") and params.get("precision"):
        props["Precision"] = params["precision"]
    return platform, props


def make_simulation(prmtop: app.AmberPrmtopFile, system: mm.System,
                    integrator: mm.Integrator, params: dict) -> app.Simulation:
    platform, props = platform_and_properties(params)
    return app.Simulation(prmtop.topology, system, integrator, platform, props)


def langevin(params: dict, temperature_K: float, gamma_ps: float, seed: int) -> mm.LangevinMiddleIntegrator:
    integ = mm.LangevinMiddleIntegrator(temperature_K * unit.kelvin,
                                        gamma_ps / unit.picosecond,
                                        float(params["dt_ps"]) * unit.picosecond)
    integ.setRandomNumberSeed(seed % 2**31)
    return integ

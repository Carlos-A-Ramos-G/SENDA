"""
senda.simulation.openmm_engine.stages
OpenMM drivers for 01_min, 02_heat, 03_equil and 04_NVT.

Every stage reads ``<replica>/openmm_params.json`` (written by senda-sim setup)
and writes AMBER-named outputs next to the AMBER engine's:

  <name>.rst7   AMBER ASCII restart (coords, velocities, box)
  <name>.xml    OpenMM state (exact continuation for the next stage)
  <name>.nc     AMBER NetCDF trajectory (wrapped, like iwrap=1)
  <name>.log    StateDataReporter CSV (AMBER ntpr)
  <name>.done   JSON completion marker (replaces "Total wall time:")

Runs are resumable: a finished sub-step (.done) is skipped, and an
interrupted MD run continues from its last checkpoint (<name>.chk.xml,
written every ntwr steps) at the step it was saved.
"""

from __future__ import annotations

import json
import logging
import math
import shutil
import time
from pathlib import Path

import numpy as np
import openmm as mm
from openmm import app, unit

from .io import (
    CheckpointReporter, NetCDFReporter, checkpoint_step, full_state, is_done,
    load_start, trim_log, write_done, write_rst7,
)
from .restraints import KCAL_TO_KJ
from .system import (
    add_barostat, add_positional, langevin, make_simulation, make_system,
    platform_and_properties,
)

log = logging.getLogger(__name__)

_KCAL_A_TO_KJ_NM = KCAL_TO_KJ * 10.0  # RMS gradient kcal/mol/A -> kJ/mol/nm

# Seed offsets keep stage/cycle/chunk/attempt random streams distinct.
_SEED_STAGE = {"heat": 1_000, "equil": 2_000, "prod": 3_000}


def _seed(params: dict, stage: str, sub: int = 0, attempt: int = 1) -> int:
    return (int(params["seed"]) + _SEED_STAGE[stage] + 10 * sub + 7919 * (attempt - 1)) % 2**31


def _positions_nm(path: Path) -> np.ndarray:
    path = Path(path)
    if path.suffix == ".xml":
        state = mm.XmlSerializer.deserialize(path.read_text())
        return state.getPositions(asNumpy=True).value_in_unit(unit.nanometer)
    return app.AmberInpcrdFile(str(path)).getPositions(asNumpy=True).value_in_unit(unit.nanometer)


def _has_velocities(path: Path) -> bool:
    return Path(path).suffix == ".xml" or app.AmberInpcrdFile(str(path)).velocities is not None


def _prefer_xml(stem: Path) -> Path:
    """Continue from <stem>.xml (exact) when present, else <stem>.rst7."""
    xml = Path(f"{stem}.xml")
    return xml if xml.exists() else Path(f"{stem}.rst7")


def _energies_kcal(state: mm.State) -> dict:
    return {
        "potential_kcal": state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole),
        "kinetic_kcal":   state.getKineticEnergy().value_in_unit(unit.kilocalorie_per_mole),
    }


def _box_A(state: mm.State) -> list[float]:
    return [float(np.linalg.norm(v.value_in_unit(unit.angstrom)))
            for v in state.getPeriodicBoxVectors()]


# ---------------------------------------------------------------------------
# Generic restartable MD run
# ---------------------------------------------------------------------------

def _run_md(params: dict, rd: Path, stem: Path, *, start: Path, nstlim: int,
            temperature: float, gamma: float, npt: bool, restraint_wt: float,
            ntpr: int, ntwx: int, ntwr: int, seed_key: tuple,
            max_attempts: int = 1, ramp: dict | None = None) -> None:
    """One MD run (one AMBER pmemd invocation), with checkpoint resume and retries.

    restraint_wt > 0 adds positional restraints referenced to *start* (the
    reference never moves on a retry -- unlike AMBER's heat_restart.in path).
    ramp = {"t_start", "ramp_end", "interval"}: linear 1 K -> T ramp over the
    first ramp_end steps, then an explicit hold at T (AMBER's &wt TEMP0).
    """
    if is_done(stem):
        log.info("%s already done -- skipping", stem.name)
        return
    chk = Path(f"{stem}.chk.xml")
    ref_nm = _positions_nm(start) if restraint_wt > 0 else None
    t0 = time.time()

    def temp_at(step: int) -> float:
        if not ramp:
            return temperature
        frac = min(step / ramp["ramp_end"], 1.0) if ramp["ramp_end"] else 1.0
        return ramp["t_start"] + (temperature - ramp["t_start"]) * frac

    for attempt in range(1, max_attempts + 1):
        resuming = chk.exists()
        step0 = checkpoint_step(chk) if resuming else 0
        seed = _seed(params, *seed_key, attempt=attempt)
        prmtop, system = make_system(params, rd, cutoff_A=float(params["cut"]))
        if ref_nm is not None:
            add_positional(system, params, rd, ref_nm, restraint_wt)
        if npt:
            add_barostat(system, {**params, "seed": seed}, temp_at(step0))
        integ = langevin(params, temp_at(step0), gamma, seed)
        sim = make_simulation(prmtop, system, integ, params)

        if resuming:
            log.info("%s: resuming from checkpoint at step %d (attempt %d)", stem.name, step0, attempt)
            load_start(sim.context, chk)
        else:
            load_start(sim.context, start)
            sim.context.applyConstraints(1e-6)
            if not _has_velocities(start):
                sim.context.setVelocitiesToTemperature(temp_at(0) * unit.kelvin, seed)
        sim.currentStep = step0

        logf = Path(f"{stem}.log")
        if resuming:
            trim_log(logf, step0)
        nc = NetCDFReporter(Path(f"{stem}.nc"), ntwx, system.getNumParticles(),
                            append=resuming, keep_frames=step0 // ntwx)
        sim.reporters += [
            app.StateDataReporter(str(logf), ntpr, step=True, time=True,
                                  potentialEnergy=True, kineticEnergy=True,
                                  totalEnergy=True, temperature=True, volume=True,
                                  density=True, speed=True,
                                  append=resuming and logf.exists()),
            nc,
            CheckpointReporter(chk, ntwr),
        ]

        try:
            seg_len = int(ramp["interval"]) if ramp else nstlim
            while sim.currentStep < nstlim:
                if ramp:
                    t = temp_at(sim.currentStep)
                    integ.setTemperature(t * unit.kelvin)
                    if npt:
                        sim.context.setParameter(mm.MonteCarloBarostat.Temperature(), t)
                sim.step(min(seg_len - sim.currentStep % seg_len, nstlim - sim.currentStep))
            state = full_state(sim.context)
            if not math.isfinite(state.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)):
                raise mm.OpenMMException("potential energy is NaN at the end of the run")
        except mm.OpenMMException as exc:
            nc.close()
            log.warning("%s: attempt %d/%d failed: %s", stem.name, attempt, max_attempts, exc)
            if attempt == max_attempts:
                raise
            continue

        nc.close()
        write_rst7(Path(f"{stem}.rst7"), state)
        Path(f"{stem}.xml").write_text(mm.XmlSerializer.serialize(state))
        write_done(stem, steps=nstlim, attempts=attempt,
                   wall_time_s=round(time.time() - t0, 1),
                   temperature_K=temperature, box_A=_box_A(state),
                   seed=seed, **_energies_kcal(state))
        for p in (chk, Path(f"{chk}.step")):
            p.unlink(missing_ok=True)
        log.info("%s done (%d steps)", stem.name, nstlim)
        return


# ---------------------------------------------------------------------------
# 01_min
# ---------------------------------------------------------------------------

def minimize(params: dict, rd: Path) -> None:
    """Repeated L-BFGS minimisations until the RMS gradient < convergence_threshold.

    Mirrors the AMBER loop (structure_min_<i>.rst7 per cycle, max_cycles_cap,
    warning fallback). No positional restraints and no constraints (AMBER
    min.in has ntr=0, ntc=1); NMR restraints are included when present.
    """
    d = rd / "01_min"
    final = d / "min_final.json"
    if final.exists():
        log.info("01_min already done -- skipping")
        return
    cfg = params["min"]
    thr_kcal = float(cfg["convergence_threshold"])
    thr_kj = thr_kcal * _KCAL_A_TO_KJ_NM

    prmtop, system = make_system(params, rd, cutoff_A=float(cfg["cut"]), constrained=False)
    platform, props = platform_and_properties(params)
    ctx = mm.Context(system, mm.VerletIntegrator(0.001), platform, props)

    start = d / "structure_min_0.rst7"
    if not start.exists():
        shutil.copy(rd / "00_prep" / "structure.rst7", start)
    done_cycles = [i for i in range(1, int(cfg["max_cycles_cap"]) + 1)
                   if is_done(d / f"structure_min_{i}")]
    last = max(done_cycles, default=0)
    load_start(ctx, d / f"structure_min_{last}.rst7")

    def _done_rms(i: int) -> float:
        return json.loads(Path(f"{d / f'structure_min_{i}'}.done").read_text())["rms_kcal_A"]

    converged = bool(last) and _done_rms(last) < thr_kcal
    first = int(cfg["max_cycles_cap"]) + 1 if converged else last + 1
    for i in range(first, int(cfg["max_cycles_cap"]) + 1):
        stem = d / f"structure_min_{i}"
        logf = Path(f"{stem}.log")
        logf.write_text(f"# cycle {i}: L-BFGS, maxIterations={cfg['maxcyc']}, "
                        f"tolerance={thr_kj:.4g} kJ/mol/nm\n")
        t0 = time.time()
        mm.LocalEnergyMinimizer.minimize(ctx, thr_kj, int(cfg["maxcyc"]))
        state = ctx.getState(getPositions=True, getEnergy=True, getForces=True,
                             enforcePeriodicBox=True)
        forces = state.getForces(asNumpy=True).value_in_unit(unit.kilojoule_per_mole / unit.nanometer)
        rms = float(np.sqrt(np.mean(forces ** 2))) / _KCAL_A_TO_KJ_NM
        epot = state.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        with logf.open("a") as fh:
            fh.write(f"ENERGY_kcal={epot:.4f} RMS_kcal_A={rms:.6g}\n")
        write_rst7(Path(f"{stem}.rst7"), state, velocities=False)
        write_done(stem, potential_kcal=epot, rms_kcal_A=rms,
                   wall_time_s=round(time.time() - t0, 1))
        log.info("01_min cycle %d: E=%.2f kcal/mol, RMS=%.3g kcal/mol/A", i, epot, rms)
        last = i
        if rms < thr_kcal:
            converged = True
            break

    if last == 0:
        raise RuntimeError("01_min: max_cycles_cap must be >= 1")
    if not converged:
        log.warning("01_min did not converge in %d cycles (RMS %.3g > %.3g); "
                    "using the last cycle", last, _done_rms(last), thr_kcal)
    final.write_text(json.dumps({"final_cycle": last, "converged": converged}, indent=2) + "\n")


# ---------------------------------------------------------------------------
# 02_heat
# ---------------------------------------------------------------------------

def heat(params: dict, rd: Path) -> None:
    final = json.loads((rd / "01_min" / "min_final.json").read_text())["final_cycle"]
    cfg = params["heat"]
    _run_md(
        params, rd, rd / "02_heat" / "structure_heat",
        start=rd / "01_min" / f"structure_min_{final}.rst7",
        nstlim=int(cfg["nstlim"]), temperature=float(params["temperature"]),
        gamma=float(cfg["gamma_ln"]), npt=True,
        restraint_wt=float(cfg["restraint_wt"]),
        ntpr=int(cfg["ntpr"]), ntwx=int(cfg["ntwx"]), ntwr=int(cfg["ntwr"]),
        seed_key=("heat",), max_attempts=int(cfg["max_retries"]),
        ramp={"t_start": float(cfg["t_start"]), "ramp_end": int(cfg["ramp_end"]),
              "interval": int(cfg["update_interval"])},
    )


# ---------------------------------------------------------------------------
# 03_equil
# ---------------------------------------------------------------------------

def equilibrate(params: dict, rd: Path) -> None:
    """Restrained NPT cycles (reference = start of each cycle), then one NVT cycle."""
    d = rd / "03_equil"
    cfg = params["equil"]
    zero = d / "structure_equil_0.rst7"
    if not zero.exists():
        shutil.copy(rd / "02_heat" / "structure_heat.rst7", zero)
    prev = _prefer_xml(rd / "02_heat" / "structure_heat")

    cycles = cfg["cycles"]
    for i, cyc in enumerate(cycles, start=1):
        stem = d / f"structure_equil_{i}"
        log.info("03_equil: cycle %d/%d (%s, wt=%s)", i, len(cycles),
                 "NPT" if cyc["npt"] else "NVT", cyc["restraint_wt"])
        _run_md(
            params, rd, stem, start=prev, nstlim=int(cyc["nstlim"]),
            temperature=float(params["temperature"]), gamma=float(cfg["gamma_ln"]),
            npt=bool(cyc["npt"]), restraint_wt=float(cyc["restraint_wt"]),
            ntpr=int(cfg["ntpr"]), ntwx=int(cfg["ntwx"]),
            ntwr=int(cfg.get("ntwr") or cyc["nstlim"]),
            seed_key=("equil", i),
        )
        prev = _prefer_xml(stem)

    last = d / f"structure_equil_{len(cycles)}"
    nvt0 = rd / "04_NVT" / "structure_NVT_0"
    for ext in (".rst7", ".xml"):
        shutil.copy(f"{last}{ext}", f"{nvt0}{ext}")


# ---------------------------------------------------------------------------
# 04_NVT
# ---------------------------------------------------------------------------

def produce(params: dict, rd: Path, chunk: int) -> None:
    cfg = params["prod"]
    d = rd / "04_NVT"
    _run_md(
        params, rd, d / f"structure_NVT_{chunk}",
        start=_prefer_xml(d / f"structure_NVT_{chunk - 1}"),
        nstlim=int(cfg["nstlim"]), temperature=float(params["temperature"]),
        gamma=float(cfg["gamma_ln"]), npt=False, restraint_wt=0.0,
        ntpr=int(cfg["ntpr"]), ntwx=int(cfg["ntwx"]), ntwr=int(cfg["ntwr"]),
        seed_key=("prod", chunk),
    )

#!/usr/bin/env python3
"""
Side-by-side validation of the AMBER and OpenMM engines on one system.

  prepare         Build two identical senda projects (engine: amber / openmm)
                  from one dimer PDB, with a shortened protocol, and run
                  `senda-sim setup` in each. Needs pmemd.cuda for the AMBER side.

  prepare-openmm  Build an OpenMM replica from the 00_prep of an AMBER run that
                  already exists (same parm7/rst7/restr, so both engines start
                  from identical inputs), using the protocol of the config.yaml
                  that run was set up with. Writes run_openmm.sh.

  energy-check    Re-evaluate frames of an existing AMBER trajectory with
                  OpenMM and with sander (AmberTools, CPU, no licence) and
                  compare the energies term by term. Validates the force field
                  on real conformations of the system, no MD needed.

  compare  Read both finished replicas and report, per engine:
             topology (atoms, residues, box)        -- must be identical
             minimisation (final energy, RMS)
             heating / NPT: temperature, density
             production: T, potential energy, CA RMSD vs the minimised
             structure, radius of gyration, inter-chain CA distance, per-residue
             CA RMSF (correlation between engines), NMR restraint values, speed
           Writes summary.md and comparison.png to --out.

Usage
-----
  python scripts/compare_engines.py prepare \\
      --pdb /path/to/protein/WT_APO_dimer.pdb \\
      --base-config /path/to/config.yaml \\
      --out engine_test --mode cluster --openmm-env "conda activate openmm"

  # submit both (cluster):  (cd engine_test/amber  && senda-sim --config config.yaml submit)
  #                         (cd engine_test/openmm && senda-sim --config config.yaml submit)

  python scripts/compare_engines.py compare --root engine_test --out engine_test/report

  # against an AMBER run that already exists (no pmemd.cuda needed):
  python scripts/compare_engines.py prepare-openmm --prep amber_run/00_prep \
      --restr amber_run/restr --config amber_run_config.yaml --out openmm_run --prod-ns 20
  python scripts/compare_engines.py energy-check --prep amber_run/00_prep \
      --traj amber_run/04_NVT/structure_NVT_1.nc --out energy_check
  bash openmm_run/run_openmm.sh
  python scripts/compare_engines.py compare --amber amber_run --openmm openmm_run --prod-ns 20

Requirements: senda importable, pyyaml, numpy, scipy, matplotlib, openmm, and
AmberTools on PATH (cpptraj for prepare-openmm/compare, sander for energy-check).
"""

from __future__ import annotations

import argparse
import copy
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import yaml

ENGINES = ("amber", "openmm")
KJ_TO_KCAL = 1.0 / 4.184

# Stage lengths per profile. "short" is the default real-system check
# (~25 ns per engine); "smoke" only exercises the plumbing.
PROFILES = {
    "short": {
        "min":        {"maxcyc": 30000, "ncyc": 500, "cut": 10.0,
                       "max_cycles_cap": 5, "convergence_threshold": 3.0e-3},
        "heat":       {"ps": 200, "max_retries": 5},
        "equil":      {"npt_ns": 0.5, "nvt_ns": 1.0},
        "production": {"ns_per_chunk": 10, "total_chunks": 2, "chunks_per_job": 2},
        "output":     {"ntpr": 5000, "ntwx": 5000, "ntwr": 50000},
    },
    "smoke": {
        "min":        {"maxcyc": 200, "ncyc": 50, "cut": 10.0,
                       "max_cycles_cap": 2, "convergence_threshold": 3.0e-3},
        "heat":       {"ps": 4, "max_retries": 2},
        "equil":      {"npt_ns": 0.004, "nvt_ns": 0.004},   # >= 1 AMBER frame (ntwx=1000)
        "production": {"ns_per_chunk": 0.004, "total_chunks": 1, "chunks_per_job": 1},
        "output":     {"ntpr": 100, "ntwx": 100, "ntwr": 500},
        "equil_schedule": [5.0, 2.0],
    },
}

# Settings copied from the base config's amber_simulator block (protocol
# lengths come from the profile instead).
_INHERITED_SIM_KEYS = ("use_hmr", "temperature", "leap", "restraints", "amber_home", "openmm")


# ===========================================================================
# prepare
# ===========================================================================

def _parse_pdb_name(pdb: Path) -> tuple[str, str]:
    m = re.fullmatch(r"(.+)_([^_]+)_dimer\.pdb", pdb.name)
    if not m:
        sys.exit(f"{pdb.name}: expected <mutant>_<inhibitor>_dimer.pdb (senda naming)")
    return m.group(1), m.group(2)


def _engine_config(base: dict, engine: str, inh: str, mut: str, args) -> dict:
    prof = PROFILES[args.profile]
    base_sim = base.get("amber_simulator") or {}
    sim = {k: copy.deepcopy(base_sim[k]) for k in _INHERITED_SIM_KEYS if k in base_sim}
    sim.setdefault("use_hmr", True)
    sim.setdefault("temperature", 300.0)
    for key in ("min", "heat", "equil", "production", "output"):
        sim[key] = dict(prof[key])
    if "equil_schedule" in prof:
        pos = sim.setdefault("restraints", {}).setdefault("positional", {}) or {}
        pos["equil_schedule"] = list(prof["equil_schedule"])
        sim["restraints"]["positional"] = pos
    sim["engine"] = engine
    sim["execution_mode"] = args.mode
    if engine == "openmm":
        omm = sim.setdefault("openmm", {}) or {}
        omm["platform"] = args.platform
        omm["seed"] = args.seed
        if args.openmm_command:
            omm["command"] = args.openmm_command
        sim["openmm"] = omm

    mc = {"output_dir": "protein"}
    disulfides = (base.get("michaelis_complex") or {}).get("disulfides")
    if disulfides:
        mc["disulfides"] = disulfides

    slurm = copy.deepcopy(base.get("slurm") or {})
    if args.openmm_env:
        target = slurm
        if slurm.get("profiles"):
            target = slurm["profiles"][slurm["active_profile"]]
        target["openmm_module"] = args.openmm_env

    return {"inhibitors": [inh], "mutants": [mut], "replicas": 1,
            "amber_simulator": sim, "michaelis_complex": mc, "slurm": slurm}


def _run_senda_sim(project: Path, *extra: str) -> None:
    cmd = [sys.executable, "-c", "from senda.simulation.run import main; main()",
           "--config", "config.yaml", *extra]
    res = subprocess.run(cmd, cwd=project, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"senda-sim {' '.join(extra)} failed in {project}:\n{res.stdout}\n{res.stderr}")


def prepare(args) -> None:
    pdb = Path(args.pdb).resolve()
    if not pdb.exists():
        sys.exit(f"PDB not found: {pdb}")
    mut, inh = _parse_pdb_name(pdb)
    base = {}
    if args.base_config:
        base = yaml.safe_load(Path(args.base_config).read_text()) or {}
    out = Path(args.out).resolve()

    lig_src = None
    if inh != "APO":
        lig_src = Path(args.ligand_lib or pdb.parent.parent / "ligands_libraries") / inh
        if not lig_src.is_dir():
            sys.exit(f"ligand parameters for {inh} not found at {lig_src} (use --ligand-lib)")

    for engine in ENGINES:
        project = out / engine
        if project.exists() and not args.force:
            sys.exit(f"{project} exists -- use --force to overwrite")
        (project / "protein").mkdir(parents=True, exist_ok=True)
        shutil.copy(pdb, project / "protein" / pdb.name)
        if lig_src:
            shutil.copytree(lig_src, project / "ligands_libraries" / inh, dirs_exist_ok=True)
        cfg = _engine_config(base, engine, inh, mut, args)
        (project / "config.yaml").write_text(
            f"# senda engine comparison ({engine}), generated by compare_engines.py\n"
            + yaml.safe_dump(cfg, sort_keys=False)
        )
        setup_args = ["setup", "--mode", args.mode, "--force"]
        if args.submit:
            setup_args.append("--submit")
        _run_senda_sim(project, *setup_args)
        print(f"  {engine:<7} {project / 'simulations' / inh / mut / 'replica_1'}")

    print(f"\nProfile '{args.profile}' written to {out}.")
    if not args.submit:
        launch = "submit" if args.mode == "cluster" else "submit --mode local"
        for engine in ENGINES:
            print(f"  (cd {out / engine} && senda-sim --config config.yaml {launch})")
    print(f"Then: python {Path(__file__).name} compare --root {out}")


# ===========================================================================
# prepare-openmm: OpenMM replica from an existing AMBER 00_prep
# ===========================================================================

RUN_OPENMM_SH = """\
#!/bin/bash
# OpenMM run of stages 01_min..04_NVT, generated by compare_engines.py prepare-openmm.
# Restartable: completed stages are skipped and an interrupted MD run resumes
# from its last checkpoint, so it is safe to resubmit after a walltime kill.
set -euo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OMM="${OMM:-__OMM__}"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

log "01_min";   $OMM run --stage min   --replica-dir "$DIR"
log "02_heat";  $OMM run --stage heat  --replica-dir "$DIR"
log "03_equil"; $OMM run --stage equil --replica-dir "$DIR"
for ((i=1; i<=__CHUNKS__; i++)); do
    log "04_NVT chunk ${i}/__CHUNKS__"
    $OMM run --stage prod --chunk "${i}" --replica-dir "$DIR"
done
log "done"
"""


def _dimer_name(prep: Path) -> tuple[str, str]:
    pdbs = sorted(prep.glob("*_dimer.pdb"))
    if not pdbs:
        sys.exit(f"no <mutant>_<inhibitor>_dimer.pdb in {prep}")
    return _parse_pdb_name(pdbs[0])


def prepare_openmm(args) -> None:
    from senda.simulation.setup import _get_nmr_list, _write_openmm_params

    prep = Path(args.prep).resolve()
    for f in ("structure.parm7", "structure.rst7"):
        if not (prep / f).exists():
            sys.exit(f"{prep / f} not found -- --prep must be the 00_prep of the AMBER run")
    mut, inh = _dimer_name(prep)
    cfg = yaml.safe_load(Path(args.config).read_text()) or {}
    sim = copy.deepcopy(cfg.get("amber_simulator") or {})
    if sim.get("use_hmr", True) and not (prep / "structure_HMR.parm7").exists():
        sys.exit(f"use_hmr is on but {prep / 'structure_HMR.parm7'} is missing")

    sim["engine"] = "openmm"
    omm = sim.get("openmm") or {}
    omm.update(platform=args.platform, seed=args.seed)
    sim["openmm"] = omm
    if args.prod_ns:
        chunk = float(args.chunk_ns or args.prod_ns)
        n = max(1, round(float(args.prod_ns) / chunk))
        sim["production"] = {**(sim.get("production") or {}), "ns_per_chunk": chunk,
                             "total_chunks": n, "chunks_per_job": n}

    rd = Path(args.out).resolve()
    if rd.exists() and any(rd.iterdir()) and not args.force:
        sys.exit(f"{rd} exists and is not empty -- use --force")
    for d in ("01_min", "02_heat", "03_equil", "04_NVT"):
        (rd / d).mkdir(parents=True, exist_ok=True)
    if (rd / "00_prep").resolve() != prep:
        shutil.copytree(prep, rd / "00_prep", dirs_exist_ok=True)

    # NMR restraints: reuse the AMBER run's restr so both engines are restrained identically
    has_nmr = bool(_get_nmr_list(sim, inh)) or bool((sim.get("restraints") or {}).get("custom_restrainer"))
    restr = Path(args.restr) if args.restr else None
    if restr:
        shutil.copy(restr, rd / "restr")
    elif has_nmr:
        sys.exit(f"the config has NMR restraints for {inh}: pass --restr <AMBER replica>/restr "
                 "so OpenMM uses exactly the restraints of the AMBER run")

    _write_openmm_params(rd, inh, 1, sim)
    params_path = rd / "openmm_params.json"
    params = json.loads(params_path.read_text())
    params["nmr"] = "restr" if (rd / "restr").exists() else None
    params_path.write_text(json.dumps(params, indent=2) + "\n")

    # positional-restraint mask -> atom list (same ambermask parser pmemd uses)
    res = subprocess.run(["cpptraj", "-i", "restraint_mask.cpptraj"], cwd=rd / "00_prep",
                         capture_output=True, text=True)
    (rd / "00_prep" / "restraint_mask.log").write_text(res.stdout + res.stderr)
    if res.returncode != 0 or not (rd / "00_prep" / "restraint_atoms.dat").exists():
        sys.exit("cpptraj could not resolve the restraint mask -- is AmberTools on PATH? "
                 f"see {rd / '00_prep' / 'restraint_mask.log'}")

    omm_cmd = args.openmm_command or omm.get("command") or "senda-omm"
    script = rd / "run_openmm.sh"
    script.write_text(RUN_OPENMM_SH.replace("__OMM__", omm_cmd)
                      .replace("__CHUNKS__", str(params["prod"]["total_chunks"])))
    script.chmod(0o755)

    n_atoms = sum(1 for _ in open(rd / "00_prep" / "restraint_atoms.dat")) - 1
    print(f"OpenMM replica: {rd}")
    print(f"  system        : {mut}_{inh}  (topology {params['topology']}, dt {params['dt_ps']} ps)")
    print(f"  NMR restraints: {'restr (from the AMBER run)' if params['nmr'] else 'none'}")
    print(f"  positional    : {params['positional_mask']!r} -> {n_atoms} atoms")
    print(f"  equilibration : {len(params['equil']['cycles'])} cycles")
    print(f"  production    : {params['prod']['total_chunks']} x "
          f"{params['prod']['nstlim'] * params['dt_ps'] / 1000:g} ns")
    print(f"  platform      : {params['platform']} ({params['precision']}), seed {params['seed']}")
    print(f"Run: bash {script}")


# ===========================================================================
# energy-check: AMBER trajectory frames re-evaluated by sander and OpenMM
# ===========================================================================

# AMBER's effective Coulomb constant (charges stored x 18.2223) vs OpenMM's CODATA value.
_COULOMB_AMBER, _COULOMB_OPENMM = 18.2223 ** 2, 332.0637133
_SANDER_KEYS = {r"\bBOND": "bond", r"\bANGLE": "angle", r"\bDIHED": "dihedral",
                r"\bVDWAALS": "vdw", r"(?<!1-4 )\bEEL": "eel", r"1-4 VDW": "vdw14",
                r"1-4 EEL": "eel14"}


def _extract_frames(trajs: list[Path], n: int, out_nc: Path, skip_ns: float) -> list[float]:
    """Write ~n evenly spaced frames of *trajs* to *out_nc*; return their times (ps)."""
    from scipy.io import netcdf_file
    from senda.simulation.openmm_engine.io import AmberNetCDFWriter

    index = []                      # (file, frame, time)
    for t in trajs:
        with netcdf_file(str(t), "r", mmap=False, maskandscale=False) as nc:
            times = np.array(nc.variables["time"].data, dtype=float)
        index += [(t, i, tm) for i, tm in enumerate(times)]
    if not index:
        sys.exit("the trajectories contain no frames")
    t0 = index[0][2]
    index = [x for x in index if x[2] - t0 >= skip_ns * 1000.0] or index
    pick = [index[int(round(k))] for k in np.linspace(0, len(index) - 1, min(n, len(index)))]

    writer = None
    for t in dict.fromkeys(p for p, _, _ in pick):
        with netcdf_file(str(t), "r", mmap=True, maskandscale=False) as nc:
            crd, lengths, angles = (nc.variables[k] for k in ("coordinates", "cell_lengths", "cell_angles"))
            for p, i, tm in pick:
                if p != t:
                    continue
                if writer is None:
                    writer = AmberNetCDFWriter(out_nc, crd.shape[1])
                writer.write_frame(np.array(crd[i]), tm, list(map(float, lengths[i])),
                                   list(map(float, angles[i])))
            del crd, lengths, angles
    writer.close()
    return [tm for _, _, tm in pick]


def _sander_energies(prep: Path, frames: Path, work: Path, cut: float, dsum_tol: float) -> list[dict]:
    (work / "energy.in").write_text(
        "energy re-evaluation of trajectory frames\n &cntrl\n"
        f"  imin=5, maxcyc=1, ntb=1, cut={cut}, ntpr=1, ntx=1,\n /\n"
        f" &ewald dsum_tol={dsum_tol}, /\n")
    res = subprocess.run(["sander", "-O", "-i", "energy.in", "-o", "energy.out",
                          "-p", str(prep / "structure.parm7"), "-c", str(prep / "structure.rst7"),
                          "-y", str(frames), "-r", "energy.rst7", "-inf", "energy.mdinfo"],
                         cwd=work, capture_output=True, text=True)
    if res.returncode != 0:
        sys.exit(f"sander failed:\n{res.stderr[-2000:]}\n(see {work / 'energy.out'})")
    out, cur = [], None
    for line in (work / "energy.out").read_text(errors="replace").splitlines():
        if "NSTEP" in line and "ENERGY" in line and "RMS" in line:
            cur = {}
            out.append(cur)
            continue
        if cur is None:
            continue
        for pattern, name in _SANDER_KEYS.items():
            m = re.search(pattern + r"\s*=\s*(-?[\d.]+|\*+)", line)
            if m:
                cur[name] = float(m.group(1)) if "*" not in m.group(1) else float("nan")
    return [r for r in out if len(r) == len(_SANDER_KEYS)]


def _openmm_energies(prep: Path, frames: Path, cut: float, ewald_tol: float, platform: str) -> list[dict]:
    import openmm as mm
    from openmm import app, unit
    from scipy.io import netcdf_file

    prmtop = app.AmberPrmtopFile(str(prep / "structure.parm7"))

    def build(scale_coulomb: bool):
        system = prmtop.createSystem(nonbondedMethod=app.PME, nonbondedCutoff=cut * 0.1 * unit.nanometer,
                                     constraints=None, rigidWater=False, ewaldErrorTolerance=ewald_tol)
        groups = {}
        for f in system.getForces():
            name = {"HarmonicBondForce": "bond", "HarmonicAngleForce": "angle",
                    "PeriodicTorsionForce": "dihedral", "NonbondedForce": "nonbonded"}.get(type(f).__name__)
            if name is None:
                continue
            f.setForceGroup(len(groups) + 1)
            groups[name] = len(groups) + 1
            if name == "nonbonded" and scale_coulomb:
                s = math.sqrt(_COULOMB_AMBER / _COULOMB_OPENMM)
                for i in range(f.getNumParticles()):
                    q, sig, eps = f.getParticleParameters(i)
                    f.setParticleParameters(i, q * s, sig, eps)
                for i in range(f.getNumExceptions()):
                    a, b, qq, sig, eps = f.getExceptionParameters(i)
                    f.setExceptionParameters(i, a, b, qq * s * s, sig, eps)
        ctx = mm.Context(system, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName(platform))
        return ctx, groups

    contexts = {"raw": build(False), "corrected": build(True)}
    out = []
    with netcdf_file(str(frames), "r", mmap=False, maskandscale=False) as nc:
        crd = np.array(nc.variables["coordinates"].data, dtype=float)
        lengths = np.array(nc.variables["cell_lengths"].data, dtype=float)
    for k in range(len(crd)):
        rec = {}
        for tag, (ctx, groups) in contexts.items():
            a, b, c = lengths[k] * 0.1
            ctx.setPeriodicBoxVectors(mm.Vec3(a, 0, 0), mm.Vec3(0, b, 0), mm.Vec3(0, 0, c))
            ctx.setPositions(crd[k] * 0.1)
            for name, g in groups.items():
                if tag == "corrected" and name != "nonbonded":
                    continue
                e = ctx.getState(getEnergy=True, groups={g}).getPotentialEnergy()
                rec[name if tag == "raw" else "nonbonded_coulomb_corrected"] = \
                    e.value_in_unit(unit.kilocalorie_per_mole)
        out.append(rec)
    return out


def energy_check(args) -> None:
    prep = Path(args.prep).resolve()
    trajs = [Path(t).resolve() for t in args.traj]
    work = Path(args.out).resolve()
    work.mkdir(parents=True, exist_ok=True)
    frames = work / "frames.nc"

    times = _extract_frames(trajs, args.nframes, frames, args.skip_ns)
    print(f"{len(times)} frames from {len(trajs)} trajectory file(s) -> {frames}")
    print("sander ...")
    amb = _sander_energies(prep, frames, work, args.cut, args.dsum_tol)
    print(f"OpenMM ({args.platform}) ...")
    omm = _openmm_energies(prep, frames, args.cut, args.ewald_tol, args.platform)
    if len(amb) != len(omm):
        sys.exit(f"sander returned {len(amb)} frames, OpenMM {len(omm)} -- see {work / 'energy.out'}")

    terms = [("Bond", "bond", "bond"), ("Angle", "angle", "angle"),
             ("Dihedral (incl. impropers)", "dihedral", "dihedral"),
             ("Nonbonded (vdW + elec + 1-4)", ("vdw", "eel", "vdw14", "eel14"), "nonbonded"),
             ("Nonbonded, Coulomb constant matched ²", ("vdw", "eel", "vdw14", "eel14"),
              "nonbonded_coulomb_corrected")]
    rows, per_frame = [], {"time_ps": times}
    for label, akey, okey in terms:
        a = np.array([sum(r[k] for k in akey) if isinstance(akey, tuple) else r[akey] for r in amb])
        o = np.array([r[okey] for r in omm])
        d = o - a
        per_frame[f"{okey}_amber"], per_frame[f"{okey}_openmm"] = a, o
        rows.append((label, a.mean(), o.mean(), d.mean(), np.abs(d).max(),
                     np.abs(d).max() / max(np.abs(a).mean(), 1e-12),
                     float(np.corrcoef(a, o)[0, 1]) if a.std() > 0 and o.std() > 0 else float("nan")))
    eel = np.array([r["eel"] + r["eel14"] for r in amb])

    with open(work / "energy_check.csv", "w") as fh:
        keys = list(per_frame)
        fh.write(",".join(keys) + "\n")
        for k in range(len(times)):
            fh.write(",".join(f"{per_frame[c][k]:.6f}" for c in keys) + "\n")

    md = ["# Energy re-evaluation: sander vs OpenMM", "",
          f"- Topology: `{prep / 'structure.parm7'}`",
          f"- Frames: {len(times)} from {', '.join(f'`{t.name}`' for t in trajs)} "
          f"({times[0] / 1000:.2f}–{times[-1] / 1000:.2f} ns)",
          f"- Cutoff {args.cut} Å, PME: sander dsum_tol {args.dsum_tol}, "
          f"OpenMM ewaldErrorTolerance {args.ewald_tol}; no restraints, no constraints", "",
          "| Term | sander mean | OpenMM mean | mean Δ | max \\|Δ\\| | max \\|Δ\\| / \\|E\\| | r (frames) |",
          "|---|---|---|---|---|---|---|"]
    md += [f"| {l} | {a:.3f} | {o:.3f} | {d:+.4f} | {mx:.4f} | {rel:.1e} | {r:.6f} |"
           for l, a, o, d, mx, rel, r in rows]
    expected = (_COULOMB_OPENMM / _COULOMB_AMBER - 1) * eel.mean()
    md += ["", "All energies in kcal/mol; Δ = OpenMM − sander.", "",
           "**How to read it**",
           "- Bond, angle and dihedral should agree to ~1e-3 kcal/mol (same parm7, same functional forms).",
           f"- Raw nonbonded differs by the Coulomb-constant convention: expected ≈ "
           f"{expected:+.2f} kcal/mol here (3.5e-5 × the electrostatic energy).",
           "- ² With OpenMM's charges rescaled to AMBER's constant, the remaining nonbonded "
           "difference is PME/cutoff numerics only: should be well below 1 kcal/mol and "
           "r ≈ 1 across frames.", ""]
    text = "\n".join(md)
    (work / "energy_check.md").write_text(text)
    print("\n" + text)
    print(f"Per-frame values: {work / 'energy_check.csv'}")


# ===========================================================================
# compare: readers
# ===========================================================================

def _natural(paths) -> list[Path]:
    def key(p):
        m = re.search(r"_(\d+)\.\w+$", p.name)
        return int(m.group(1)) if m else -1
    return sorted(paths, key=key)


def _replica(project: Path) -> Path:
    """A replica directory (contains 00_prep) or a senda project holding one."""
    if (project / "00_prep").is_dir():
        return project
    found = sorted(project.glob("simulations/*/*/replica_1"))
    if not found:
        sys.exit(f"no simulations/*/*/replica_1 under {project}")
    return found[0]


_RECORD_KV = re.compile(r"([A-Za-z0-9\-\(\)][A-Za-z0-9\-\(\) ]*?)\s*=\s*(-?\d[\d.]*(?:[eE][-+]?\d+)?)")


def read_mdout(path: Path) -> list[dict]:
    """Per-step records of an AMBER MD .out (stops at the AVERAGES block)."""
    records, cur = [], None
    in_results = False
    for line in path.read_text(errors="replace").splitlines():
        if "4.  RESULTS" in line:
            in_results = True
            continue
        if not in_results:
            continue
        if "A V E R A G E S" in line:
            break
        if line.lstrip().startswith("NSTEP ="):
            cur = {}
            records.append(cur)
        if cur is not None:
            for k, v in _RECORD_KV.findall(line):
                cur[k.strip()] = float(v)
    out = []
    for r in records:
        if "TIME(PS)" not in r:
            continue
        out.append({"time": r["TIME(PS)"], "T": r.get("TEMP(K)"), "epot": r.get("EPtot"),
                    "etot": r.get("Etot"), "density": r.get("Density"),
                    "volume": r.get("VOLUME") or None})
    return out


def read_omm_log(path: Path) -> list[dict]:
    """Records of an OpenMM StateDataReporter CSV, in AMBER units (kcal/mol, A^3)."""
    lines = path.read_text().splitlines()
    if not lines:
        return []
    cols = [c.strip('#" ') for c in lines[0].split(",")]
    idx = {name: i for i, name in enumerate(cols)}

    def col(row, name, scale=1.0):
        i = idx.get(name)
        return float(row[i]) * scale if i is not None else None

    out = []
    for line in lines[1:]:
        row = line.split(",")
        try:
            out.append({
                "time": col(row, "Time (ps)"), "T": col(row, "Temperature (K)"),
                "epot": col(row, "Potential Energy (kJ/mole)", KJ_TO_KCAL),
                "etot": col(row, "Total Energy (kJ/mole)", KJ_TO_KCAL),
                "density": col(row, "Density (g/mL)"),
                "volume": col(row, "Box Volume (nm^3)", 1000.0),
                "speed": col(row, "Speed (ns/day)"),
            })
        except (ValueError, IndexError):
            continue
    return out


def _amber_speed(path: Path) -> float | None:
    vals = re.findall(r"ns/day =\s*([\d.]+)", path.read_text(errors="replace"))
    return float(vals[-1]) if vals else None


def _amber_min(path: Path) -> dict:
    lines = path.read_text(errors="replace").splitlines()
    for i, line in enumerate(lines):
        if "FINAL RESULTS" in line:
            for j in range(i, min(i + 10, len(lines))):
                if "NSTEP" in lines[j] and "ENERGY" in lines[j]:
                    f = lines[j + 1].split()
                    return {"energy_kcal": float(f[1]), "rms_kcal_A": float(f[2])}
    return {}


def stage_runs(rd: Path, engine: str) -> list[tuple[str, Path]]:
    """(label, per-run log) in simulation order for heat, equil and production."""
    ext = ".out" if engine == "amber" else ".log"
    runs = [("heat", rd / "02_heat" / f"structure_heat{ext}")]
    runs += [(f"equil_{p.stem.rsplit('_', 1)[1]}", p)
             for p in _natural((rd / "03_equil").glob(f"structure_equil_[1-9]*{ext}"))]
    runs += [(f"prod_{p.stem.rsplit('_', 1)[1]}", p)
             for p in _natural((rd / "04_NVT").glob(f"structure_NVT_[1-9]*{ext}"))]
    return [(lab, p) for lab, p in runs if p.exists()]


def _openmm_nvt_labels(rd: Path) -> set[str]:
    """Stages run at constant volume (OpenMM logs a constant density there; AMBER logs none)."""
    params = json.loads((rd / "openmm_params.json").read_text())
    labels = {f"equil_{i}" for i, c in enumerate(params["equil"]["cycles"], start=1) if not c["npt"]}
    return labels | {f"prod_{i}" for i in range(1, int(params["prod"]["total_chunks"]) + 1)}


def read_thermo(rd: Path, engine: str) -> dict:
    series, speeds = {}, {}
    nvt = _openmm_nvt_labels(rd) if engine == "openmm" else set()
    for label, path in stage_runs(rd, engine):
        recs = read_mdout(path) if engine == "amber" else read_omm_log(path)
        if label in nvt:
            for r in recs:
                r["density"] = None
        series[label] = recs
        if engine == "amber":
            speeds[label] = _amber_speed(path)
        else:
            sp = [r["speed"] for r in recs if r.get("speed")]
            speeds[label] = sp[-1] if sp else None
    return {"series": series, "speed": speeds}


def read_min(rd: Path, engine: str) -> dict:
    d = rd / "01_min"
    if engine == "amber":
        outs = _natural(d.glob("structure_min_[1-9]*.out"))
        if not outs:
            return {}
        return {"cycles": len(outs), **_amber_min(outs[-1]), "final": outs[-1].stem}
    final = json.loads((d / "min_final.json").read_text())
    n = final["final_cycle"]
    done = json.loads((d / f"structure_min_{n}.done").read_text())
    return {"cycles": n, "energy_kcal": done["potential_kcal"], "rms_kcal_A": done["rms_kcal_A"],
            "converged": final["converged"], "final": f"structure_min_{n}"}


# ---------------------------------------------------------------------------
# Topology
# ---------------------------------------------------------------------------

def _prmtop_pointers(path: Path) -> dict:
    text = path.read_text()
    m = re.search(r"%FLAG POINTERS\s*\n%FORMAT\(10I8\)\s*\n(.*?)%FLAG", text, re.DOTALL)
    vals = [int(x) for x in m.group(1).split()]
    return {"atoms": vals[0], "residues": vals[11]}


def _rst7_box(path: Path) -> list[float]:
    return [float(x) for x in path.read_text().splitlines()[-1].split()[:3]]


def topology_info(rd: Path) -> dict:
    prep = rd / "00_prep"
    info = _prmtop_pointers(prep / "structure.parm7")
    info["box_A"] = _rst7_box(prep / "structure.rst7")
    return info


def chain_ca_masks(rd: Path) -> list[str]:
    """CA masks per PDB chain, from the dimer PDB in 00_prep (tleap keeps residue order)."""
    from senda.analysis.distances import _parse_pdb_residues
    pdb = next((rd / "00_prep").glob("*_dimer.pdb"))
    residues = _parse_pdb_residues(pdb)
    masks, start = [], 1
    chains = []
    for ch, _, _ in residues:
        if not chains or chains[-1][0] != ch:
            chains.append([ch, 0])
        chains[-1][1] += 1
    for ch, n in chains:
        masks.append((ch, f":{start}-{start + n - 1}@CA"))
        start += n
    return masks


def nmr_restraints(rd: Path) -> list[dict]:
    path = rd / "restr"
    if not path.exists():
        return []
    out = []
    for block in re.findall(r"&rst(.*?)/", path.read_text(), re.DOTALL):
        iat = [int(x) for x in re.search(r"iat\s*=\s*([\d,\s]+)", block).group(1).replace(",", " ").split() if int(x) > 0]
        vals = {k: float(v) for k, v in re.findall(r"(r[1-4])\s*=\s*(-?[\d.]+)", block)}
        out.append({"atoms": iat, **vals})
    return out


# ---------------------------------------------------------------------------
# Structural analysis (cpptraj)
# ---------------------------------------------------------------------------

def _nframes(path: Path) -> int:
    from scipy.io import netcdf_file
    with netcdf_file(str(path), "r", mmap=False, maskandscale=False) as nc:
        return int(nc.variables["time"].shape[0])


def _trajectories(rd: Path) -> tuple[list[Path], list[Path]]:
    """Heat+equil and production trajectories, skipping files with no frames
    (a run shorter than its ntwx writes an empty .nc, which cpptraj rejects)."""
    pre = [rd / "02_heat" / "structure_heat.nc"]
    pre += _natural((rd / "03_equil").glob("structure_equil_[1-9]*.nc"))
    prod = _natural((rd / "04_NVT").glob("structure_NVT_[1-9]*.nc"))
    keep = lambda paths: [p for p in paths if p.exists() and _nframes(p) > 0]
    return keep(pre), keep(prod)


def _nc_times(paths: list[Path]) -> np.ndarray:
    from scipy.io import netcdf_file
    times = []
    for p in paths:
        with netcdf_file(str(p), "r", mmap=False, maskandscale=False) as nc:
            times.append(np.array(nc.variables["time"].data, dtype=float))
    return np.concatenate(times) if times else np.array([])


def _read_dat(path: Path) -> np.ndarray:
    rows = [l.split() for l in path.read_text().splitlines() if l.strip() and not l.startswith("#")]
    return np.array(rows, dtype=float) if rows else np.empty((0, 2))


def _prod_window(prod: list[Path], prod_ns: float | None) -> list[tuple[Path, int]]:
    """(file, last frame, 1-based) covering the first *prod_ns* ns of production."""
    if not prod_ns:
        return [(p, _nframes(p)) for p in prod]
    t0, out = None, []
    for p in prod:
        t = _nc_times([p])
        t0 = t[0] if t0 is None else t0
        n = int(np.sum(t - t0 <= prod_ns * 1000.0 + 1e-6))
        if n:
            out.append((p, n))
        if n < len(t):
            break
    return out


def run_cpptraj(rd: Path, workdir: Path, min_final: str, restraints: list[dict],
                prod_ns: float | None = None) -> dict:
    workdir.mkdir(parents=True, exist_ok=True)
    pre, prod_all = _trajectories(rd)
    window = _prod_window(prod_all, prod_ns)
    prod = [p for p, _ in window]
    if not prod:
        sys.exit(f"no production trajectories in {rd / '04_NVT'}")
    masks = chain_ca_masks(rd)
    lines = [f"parm {rd / '00_prep' / 'structure.parm7'}",
             f"reference {rd / '01_min' / (min_final + '.rst7')} [minref]"]
    lines += [f"trajin {p}" for p in pre]
    lines += [f"trajin {p} 1 {n}" for p, n in window]
    lines += ["autoimage",
              "rms rmsd @CA ref [minref] out rmsd.dat",
              "radgyr rg @CA out rg.dat"]
    if len(masks) == 2:
        lines.append(f"distance chains {masks[0][1]} {masks[1][1]} out chaindist.dat")
    kinds = {2: "distance", 3: "angle", 4: "dihedral"}
    for i, r in enumerate(restraints):
        atoms = " ".join(f"@{a}" for a in r["atoms"])
        lines.append(f"{kinds[len(r['atoms'])]} nmr{i} {atoms} out nmr{i}.dat")
    lines += ["run", "clear trajin"]
    lines += [f"trajin {p} 1 {n}" for p, n in window]
    lines += ["autoimage", "rms first @CA", "atomicfluct rmsf @CA byres out rmsf.dat", "run", "quit"]
    (workdir / "analysis.cpptraj").write_text("\n".join(lines) + "\n")
    res = subprocess.run(["cpptraj", "-i", "analysis.cpptraj"], cwd=workdir,
                         capture_output=True, text=True)
    (workdir / "cpptraj.log").write_text(res.stdout + res.stderr)
    if res.returncode != 0:
        sys.exit(f"cpptraj failed -- see {workdir / 'cpptraj.log'}")

    t_prod = np.concatenate([_nc_times([p])[:n] for p, n in window])
    t_all = np.concatenate([_nc_times(pre), t_prod]) if pre else t_prod
    t_prod0 = t_prod[0]
    data = {"time": t_all, "prod_start": t_prod0}
    for name in ("rmsd", "rg", "chaindist"):
        f = workdir / f"{name}.dat"
        if f.exists():
            data[name] = _read_dat(f)[:, 1]
    data["nmr"] = [_read_dat(workdir / f"nmr{i}.dat")[:, 1] for i in range(len(restraints))]
    rmsf = _read_dat(workdir / "rmsf.dat")
    data["rmsf"] = rmsf[:, 1] if len(rmsf) else np.array([])
    data["rmsf_res"] = rmsf[:, 0] if len(rmsf) else np.array([])
    return data


# ===========================================================================
# compare: summary
# ===========================================================================

def _ms(x) -> str:
    x = np.asarray([v for v in x if v is not None], dtype=float)
    if not len(x):
        return "--"
    return f"{x.mean():.3f} ± {x.std():.3f}" if abs(x.mean()) < 100 else f"{x.mean():.1f} ± {x.std():.1f}"


def _concat(series: dict, prefix: str, key: str) -> list:
    return [r[key] for lab, recs in series.items() if lab.startswith(prefix) for r in recs]


def _last_npt(series: dict) -> list[dict]:
    npt = [lab for lab, recs in series.items()
           if lab.startswith("equil") and any(r.get("density") for r in recs)]
    if not npt:
        return []
    recs = series[npt[-1]]
    return recs[len(recs) // 2:]


def _speed(speeds: dict, prefix: str) -> str:
    vals = [v for k, v in speeds.items() if k.startswith(prefix) and v]
    return f"{np.mean(vals):.1f}" if vals else "--"


def build_summary(res: dict) -> list[tuple[str, str, str]]:
    a, o = res["amber"], res["openmm"]
    rows = []
    rows.append(("Atoms / residues", f"{a['topo']['atoms']} / {a['topo']['residues']}",
                 f"{o['topo']['atoms']} / {o['topo']['residues']}"))
    rows.append(("Initial box (Å)", " × ".join(f"{x:.2f}" for x in a["topo"]["box_A"]),
                 " × ".join(f"{x:.2f}" for x in o["topo"]["box_A"])))
    rows.append(("Min: cycles", str(a["min"].get("cycles", "--")),
                 f"{o['min'].get('cycles', '--')} ({'converged' if o['min'].get('converged') else 'not converged'})"))
    rows.append(("Min: final E (kcal/mol)", f"{a['min'].get('energy_kcal', float('nan')):.1f}",
                 f"{o['min'].get('energy_kcal', float('nan')):.1f}"))
    rows.append(("Min: final RMS grad (kcal/mol/Å)", f"{a['min'].get('rms_kcal_A', float('nan')):.4g}",
                 f"{o['min'].get('rms_kcal_A', float('nan')):.4g}"))
    for eng in (a, o):
        heat = eng["thermo"]["series"].get("heat", [])
        eng["_heat_tail"] = heat[int(0.8 * len(heat)):] if heat else []
    rows.append(("Heat: T last 20% (K)", _ms(r["T"] for r in a["_heat_tail"]),
                 _ms(r["T"] for r in o["_heat_tail"])))
    rows.append(("NPT: density, last cycle 2nd half (g/mL)",
                 _ms(r["density"] for r in _last_npt(a["thermo"]["series"])),
                 _ms(r["density"] for r in _last_npt(o["thermo"]["series"]))))
    rows.append(("Prod: T (K)", _ms(_concat(a["thermo"]["series"], "prod", "T")),
                 _ms(_concat(o["thermo"]["series"], "prod", "T"))))
    rows.append(("Prod: E_pot (kcal/mol) ¹", _ms(_concat(a["thermo"]["series"], "prod", "epot")),
                 _ms(_concat(o["thermo"]["series"], "prod", "epot"))))
    for key, label in (("rmsd", "Prod: CA RMSD vs minimised (Å)"), ("rg", "Prod: CA Rg (Å)"),
                       ("chaindist", "Prod: chain A–B CA COM distance (Å)")):
        if key in a["struct"] and key in o["struct"]:
            vals = []
            for eng in (a, o):
                s = eng["struct"]
                vals.append(_ms(s[key][s["time"] >= s["prod_start"]]))
            rows.append((label, *vals))
    ra, ro = a["struct"]["rmsf"], o["struct"]["rmsf"]
    if len(ra) and len(ra) == len(ro):
        r = float(np.corrcoef(ra, ro)[0, 1])
        rows.append(("Prod: CA RMSF mean (Å)", f"{ra.mean():.3f}", f"{ro.mean():.3f}"))
        rows.append(("Prod: CA RMSF Pearson r (AMBER vs OpenMM)", f"{r:.3f}", ""))
    for i, rst in enumerate(res["restraints"]):
        vals = []
        for eng in (a, o):
            s = eng["struct"]
            x = s["nmr"][i][s["time"] >= s["prod_start"]]
            lo, hi = rst.get("r2"), rst.get("r3")
            inside = np.mean((x >= lo) & (x <= hi)) * 100 if lo is not None else float("nan")
            vals.append(f"{_ms(x)}  ({inside:.0f}% in [{lo:g}, {hi:g}])")
        rows.append((f"Prod: NMR restraint {i} atoms {rst['atoms']}", *vals))
    rows.append(("Speed, production (ns/day)", _speed(a["thermo"]["speed"], "prod"),
                 _speed(o["thermo"]["speed"], "prod")))
    return rows


def write_summary(rows, out: Path, res: dict) -> str:
    same_topo = (res["amber"]["topo"]["atoms"] == res["openmm"]["topo"]["atoms"]
                 and res["amber"]["topo"]["residues"] == res["openmm"]["topo"]["residues"])
    md = ["# AMBER vs OpenMM engine comparison", "",
          f"- AMBER replica: `{res['amber']['rd']}`",
          f"- OpenMM replica: `{res['openmm']['rd']}`",
          f"- Topologies identical (atoms/residues): **{'yes' if same_topo else 'NO'}**",
          f"- Production window: {res.get('window', '').strip(' ()') or 'all frames'}", "",
          "| Quantity | AMBER | OpenMM |", "|---|---|---|"]
    md += [f"| {q} | {a} | {o} |" for q, a, o in rows]
    md += ["", "Values are mean ± standard deviation over the stated window.", "",
           "¹ Absolute potential energies are not expected to match exactly: OpenMM uses the "
           "CODATA Coulomb constant (332.0637 vs AMBER's 18.2223² = 332.0522), which shifts "
           "electrostatics by ~3.5e-5 relative (several kcal/mol for a solvated protein). "
           "Compare fluctuations and trends, not the offset.",
           "", "![comparison](comparison.png)", ""]
    text = "\n".join(md)
    (out / "summary.md").write_text(text)
    return text


# ===========================================================================
# compare: figure (one measure per panel; AMBER blue, OpenMM orange)
# ===========================================================================

COLORS = {"amber": "#2a78d6", "openmm": "#eb6834"}
LABELS = {"amber": "AMBER (pmemd)", "openmm": "OpenMM"}
INK, INK2, MUTED, GRID, AXIS, SURFACE = "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7", "#fcfcfb"


def plot(res: dict, out: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.size": 9, "axes.edgecolor": AXIS, "axes.labelcolor": INK2,
        "axes.titlecolor": INK, "axes.titlesize": 10, "axes.titleweight": "bold",
        "axes.titlelocation": "left", "xtick.color": MUTED, "ytick.color": MUTED,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "legend.frameon": False,
    })
    panels = [("T", "Temperature (K)"), ("density", "Density (g/mL)"),
              ("epot", "Potential energy (kcal/mol)"), ("rmsd", "CA RMSD vs minimised (Å)"),
              ("rg", "CA radius of gyration (Å)"), ("rmsf", "Production CA RMSF (Å)")]
    nmr = res["restraints"]
    panels += [(f"nmr{i}", f"NMR restraint {i} ({'Å' if len(r['atoms']) == 2 else '°'})")
               for i, r in enumerate(nmr)]
    ncol = 2
    nrow = math.ceil(len(panels) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(11, 3.0 * nrow), constrained_layout=True)
    axes = np.atleast_1d(axes).ravel()

    for ax, (key, title) in zip(axes, panels):
        ax.set_title(title)
        for eng in ENGINES:
            e = res[eng]
            kw = dict(color=COLORS[eng], lw=1.5, label=LABELS[eng])
            if key in ("T", "density", "epot"):
                t, y = [], []
                for recs in e["thermo"]["series"].values():
                    for r in recs:
                        if r.get(key) is not None:
                            t.append(r["time"])
                            y.append(r[key])
                ax.plot(np.array(t) / 1000.0, y, **kw)
                ax.set_xlabel("Time (ns)")
            elif key == "rmsf":
                s = e["struct"]
                ax.plot(s["rmsf_res"], s["rmsf"], **kw)
                ax.set_xlabel("Residue")
            else:
                s = e["struct"]
                y = s["nmr"][int(key[3:])] if key.startswith("nmr") else s.get(key)
                if y is not None and len(y):
                    ax.plot(s["time"][: len(y)] / 1000.0, y, **kw)
                ax.set_xlabel("Time (ns)")
        if key not in ("rmsf",):
            p0 = res["openmm"]["struct"]["prod_start"] / 1000.0
            ax.axvline(p0, color=AXIS, lw=1, ls="--")
        if key.startswith("nmr"):
            r = nmr[int(key[3:])]
            ax.axhspan(r["r2"], r["r3"], color=GRID, alpha=0.6, lw=0)
    for ax in axes[len(panels):]:
        ax.set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper right", ncol=2, fontsize=9, labelcolor=INK)
    fig.suptitle("AMBER vs OpenMM — dashed line marks the start of production",
                 x=0.01, ha="left", fontsize=11, color=INK)
    fig.savefig(out / "comparison.png", dpi=150)
    plt.close(fig)


def compare(args) -> None:
    root = Path(args.root).resolve() if args.root else None
    projects = {"amber": Path(args.amber) if args.amber else root / "amber",
                "openmm": Path(args.openmm) if args.openmm else root / "openmm"}
    out = Path(args.out or (root / "report" if root else "engine_report")).resolve()
    window = f" (first {args.prod_ns:g} ns of production)" if args.prod_ns else ""
    out.mkdir(parents=True, exist_ok=True)

    res = {}
    for eng, project in projects.items():
        rd = _replica(project.resolve())
        print(f"Reading {eng}: {rd}")
        mn = read_min(rd, eng)
        res[eng] = {"rd": rd, "topo": topology_info(rd), "min": mn,
                    "thermo": read_thermo(rd, eng)}
    if args.prod_ns:
        for eng in ENGINES:
            series = res[eng]["thermo"]["series"]
            prod = [r["time"] for lab, recs in series.items() if lab.startswith("prod") for r in recs]
            if prod:
                t_end = min(prod) + args.prod_ns * 1000.0
                for lab in [l for l in series if l.startswith("prod")]:
                    series[lab] = [r for r in series[lab] if r["time"] <= t_end + 1e-6]
    restraints = nmr_restraints(res["openmm"]["rd"]) or nmr_restraints(res["amber"]["rd"])
    res["restraints"] = restraints
    for eng in ENGINES:
        print(f"Running cpptraj for {eng} ...")
        res[eng]["struct"] = run_cpptraj(res[eng]["rd"], out / eng, res[eng]["min"]["final"],
                                         restraints, args.prod_ns)

    rows = build_summary(res)
    res["window"] = window
    text = write_summary(rows, out, res)
    plot(res, out)
    print("\n" + text)
    print(f"Report: {out / 'summary.md'}\nFigure: {out / 'comparison.png'}")


# ===========================================================================

def main(argv=None) -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("prepare", help="Build the AMBER and OpenMM test projects")
    p.add_argument("--pdb", required=True, help="<mutant>_<inhibitor>_dimer.pdb")
    p.add_argument("--base-config", help="config.yaml to inherit leap/restraints/slurm from")
    p.add_argument("--out", default="engine_test")
    p.add_argument("--profile", choices=sorted(PROFILES), default="short")
    p.add_argument("--mode", choices=["cluster", "local"], default="cluster")
    p.add_argument("--platform", default="CUDA", help="OpenMM platform (CUDA, OpenCL, CPU)")
    p.add_argument("--seed", type=int, default=2026, help="OpenMM seed (reproducible test)")
    p.add_argument("--openmm-env", help="shell line activating the OpenMM env (slurm.openmm_module)")
    p.add_argument("--openmm-command", help="override amber_simulator.openmm.command")
    p.add_argument("--ligand-lib", help="ligands_libraries dir (non-APO systems)")
    p.add_argument("--submit", action="store_true", help="submit both runs after setup")
    p.add_argument("--force", action="store_true", help="overwrite existing projects")

    po = sub.add_parser("prepare-openmm", help="OpenMM replica from an existing AMBER 00_prep")
    po.add_argument("--prep", required=True, help="00_prep directory of the AMBER run")
    po.add_argument("--config", required=True, help="config.yaml the AMBER run was set up with")
    po.add_argument("--restr", help="restr (DISANG) file from the AMBER replica root")
    po.add_argument("--out", default="openmm_run", help="OpenMM replica directory to create")
    po.add_argument("--prod-ns", type=float, help="total production (default: as in the config)")
    po.add_argument("--chunk-ns", type=float, help="production chunk length (default: --prod-ns)")
    po.add_argument("--platform", default="CUDA", help="OpenMM platform (CUDA, OpenCL, CPU)")
    po.add_argument("--seed", type=int, default=2026)
    po.add_argument("--openmm-command", help="command running senda-omm (default: senda-omm)")
    po.add_argument("--force", action="store_true")

    ec = sub.add_parser("energy-check", help="Re-evaluate AMBER frames with sander and OpenMM")
    ec.add_argument("--prep", required=True, help="00_prep directory of the AMBER run")
    ec.add_argument("--traj", required=True, nargs="+", help="AMBER .nc trajectory file(s)")
    ec.add_argument("--out", default="energy_check")
    ec.add_argument("--nframes", type=int, default=50)
    ec.add_argument("--skip-ns", type=float, default=0.0, help="skip the first N ns of the trajectories")
    ec.add_argument("--cut", type=float, default=10.0, help="cutoff in A (pmemd runs used 10)")
    ec.add_argument("--dsum-tol", type=float, default=1e-6)
    ec.add_argument("--ewald-tol", type=float, default=1e-6)
    ec.add_argument("--platform", default="CPU", help="OpenMM platform for the single points")

    c = sub.add_parser("compare", help="Compare two finished replicas")
    c.add_argument("--root", help="directory written by prepare (contains amber/ and openmm/)")
    c.add_argument("--amber", help="AMBER replica or project dir (overrides --root)")
    c.add_argument("--openmm", help="OpenMM replica or project dir (overrides --root)")
    c.add_argument("--out", help="report directory (default: <root>/report)")
    c.add_argument("--prod-ns", type=float,
                   help="compare only the first N ns of production of each engine")

    args = parser.parse_args(argv)
    if args.command == "prepare":
        prepare(args)
    elif args.command == "prepare-openmm":
        prepare_openmm(args)
    elif args.command == "energy-check":
        energy_check(args)
    else:
        if not args.root and not (args.amber and args.openmm):
            parser.error("compare needs --root or both --amber and --openmm")
        compare(args)


if __name__ == "__main__":
    main()

"""
senda.simulation.setup
Build simulation directory trees and submit/launch jobs.

Directory layout produced:
  simulations/{inhibitor}/{mutant}/replica_{n}/
    00_prep/   -PDB, FF parameters, leap_structure, HMR.ccptraj
    01_min/    -min.in
    02_heat/   -heat.in
    03_equil/  -equil_1.in to equil_6.in
    04_NVT/    -prod.in, run_NVT_*.cmd  (cluster) or run_NVT_*_local.sh  (local)
    run_gpu    -master SLURM script  (cluster mode)
    run_local  -master bash script   (local mode)

Protein PDBs are read from {protein_dir}/{mutant}_{inhibitor}_dimer.pdb.
Ligand FF parameters are read from {ligands_lib_dir}/{inhibitor}/{inhibitor}.lib
and .frcmod.  APO systems have no ligand parameters.

Modes
-----
  cluster  (default) -SLURM; master job runs stages 00-03, chains NVT chunks via sbatch.
  local    -No SLURM; generates run_local for a single-GPU workstation.

Restraints
----------
  Positional restraints (backbone) are applied during heat and NPT equilibration.
  Their mask and force constants can be tuned via amber_simulator.restraints.positional.

  NMR restraints (distance, angle, dihedral) are optional. When present in
  amber_simulator.restraints.nmr, a restrainer.py script is generated that
  builds the AMBER DISANG file at job time after tleap creates the topology.
  nmr may be a plain list (applies to every ligand) or a dict keyed by
  inhibitor name (per-ligand restraints).
"""

from __future__ import annotations

import math
import shutil
import subprocess
from pathlib import Path

from .templates import (
    MIN, HEAT, EQUIL_NPT, EQUIL_NVT, PROD,
    HMR_CCPTRAJ, make_restrainer,
    RUN_GPU_HEADER, RUN_LOCAL_HEADER,
    NVT_CLUSTER_HEADER, NVT_LOCAL_HEADER,
    RUN_BODY, NVT_JOB_BODY,
    fill,
)

_SUBDIRS = ("00_prep", "01_min", "02_heat", "03_equil", "04_NVT")

_DEFAULT_RESTRAINT_MASK     = "@CA,C,N,O,H &!:WAT"
_DEFAULT_HEATING_WT         = 20.0
_DEFAULT_EQUIL_SCHEDULE     = [15.0, 12.0, 9.0, 6.0, 3.0]


# ---------------------------------------------------------------------------
# tleap input
# ---------------------------------------------------------------------------

def _render_leap(inh: str, mut: str, has_ligand: bool, leap_cfg: dict) -> str:
    """
    Build the tleap input file.

    Defaults (ff14SB + TIP3P + GAFF2, Na+ neutralization, 12 A TIP3PBOX) can
    be overridden via the optional amber_simulator.leap section in config.yaml:

      amber_simulator:
        leap:
          forcefields: [leaprc.protein.ff14SB, leaprc.water.tip3p, leaprc.gaff2]
          ions: ["Na+ 0", "Cl- 0"]
          box_type: TIP3PBOX
          box_size: 12
    """
    default_ff = ["leaprc.protein.ff14SB", "leaprc.water.tip3p"]
    if has_ligand:
        default_ff.append("leaprc.gaff2")

    forcefields = leap_cfg.get("forcefields", default_ff)
    ions        = leap_cfg.get("ions", ["Na+ 0"])
    box_type    = leap_cfg.get("box_type", "TIP3PBOX")
    box_size    = leap_cfg.get("box_size", 12)

    lines = [f"source {ff}" for ff in forcefields]

    if has_ligand:
        lines += [
            f"loadoff {inh}.lib",
            f"loadamberparams {inh}.frcmod",
        ]

    lines.append(f"structure = loadpdb {mut}_{inh}_dimer.pdb")
    for ion in ions:
        lines.append(f"addions structure {ion}")
    lines.append(f"solvatebox structure {box_type} {box_size}")
    lines.append("saveamberparm structure structure.parm7 structure.rst7")
    lines.append("quit")
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------------
# NMR restraint resolution
# ---------------------------------------------------------------------------

def _get_nmr_list(sim: dict, inh: str, has_ligand: bool) -> list:
    """Return the NMR restraint list for *inh*, supporting both formats:
      - list  -> applies to every ligand (legacy)
      - dict  -> keyed by inhibitor name (per-ligand)
    """
    if not has_ligand:
        return []
    raw = (sim.get("restraints") or {}).get("nmr")
    if not raw:
        return []
    if isinstance(raw, dict):
        return raw.get(inh) or []
    return raw  # plain list: same restraints for all ligands


# ---------------------------------------------------------------------------
# AMBER input files
# ---------------------------------------------------------------------------

def _write_input_files(replica_dir: Path, inh: str, sim: dict,
                       has_ligand: bool = True) -> str:
    """Write all AMBER .in files; return the topology filename."""
    use_hmr  = bool(sim.get("use_hmr", True))
    dt       = 0.004 if use_hmr else 0.002
    topology = "structure_HMR.parm7" if use_hmr else "structure.parm7"
    temp     = float(sim["temperature"])

    restr_cfg  = sim.get("restraints", {}) or {}
    nmr_list   = _get_nmr_list(sim, inh, has_ligand)
    has_nmr    = bool(nmr_list)

    pos_cfg     = restr_cfg.get("positional", {}) or {}
    restr_mask  = pos_cfg.get("mask",            _DEFAULT_RESTRAINT_MASK)
    heating_wt  = pos_cfg.get("heating_weight",  _DEFAULT_HEATING_WT)
    equil_sched = pos_cfg.get("equil_schedule",  _DEFAULT_EQUIL_SCHEDULE)

    nmropt   = "\n  nmropt = 1," if has_nmr else ""
    wtdisang = "&wt type='END'/\nDISANG = ../restr\n" if has_nmr else ""
    disang   = "DISANG = ../restr" if has_nmr else ""

    min_cfg = sim["min"]
    (replica_dir / "01_min" / "min.in").write_text(fill(
        MIN,
        MAXCYC=min_cfg["maxcyc"], NCYC=min_cfg["ncyc"], CUT=min_cfg["cut"],
        NMROPT=nmropt, WTBLOCK=wtdisang,
    ))

    heat_steps = int(round(float(sim["heat"]["ps"]) / dt))
    ramp_end   = int(round(0.8 * heat_steps))
    (replica_dir / "02_heat" / "heat.in").write_text(fill(
        HEAT,
        TEMP=temp, NSTLIM=heat_steps, DT=dt, RAMPEND=ramp_end,
        RESTRAINTMASK=restr_mask, RESTRAINT_WT=heating_wt,
        DISANG=disang,
    ))

    npt_steps = int(round(float(sim["equil"]["npt_ns"]) * 1000.0 / dt))
    for cycle, wt in enumerate(equil_sched, start=1):
        (replica_dir / "03_equil" / f"equil_{cycle}.in").write_text(fill(
            EQUIL_NPT,
            CYCLE=cycle, NSTLIM=npt_steps, DT=dt, TEMP=temp,
            RESTRAINTMASK=restr_mask, WT=f"{wt:.1f}",
            NMROPT=nmropt, WTDISANG=wtdisang,
        ))

    nvt_steps = int(round(float(sim["equil"]["nvt_ns"]) * 1000.0 / dt))
    (replica_dir / "03_equil" / "equil_6.in").write_text(fill(
        EQUIL_NVT,
        NSTLIM=nvt_steps, DT=dt, TEMP=temp,
        NMROPT=nmropt, WTDISANG=wtdisang,
    ))

    prod       = sim["production"]
    out_cfg    = sim["output"]
    prod_steps = int(round(float(prod["ns_per_chunk"]) * 1000.0 / dt))
    (replica_dir / "04_NVT" / "prod.in").write_text(fill(
        PROD,
        NSTLIM=prod_steps, DT=dt, TEMP=temp,
        NTPR=out_cfg["ntpr"], NTWX=out_cfg["ntwx"], NTWR=out_cfg["ntwr"],
        NMROPT=nmropt, WTDISANG=wtdisang,
    ))

    return topology


# ---------------------------------------------------------------------------
# Run scripts (cluster and local)
# ---------------------------------------------------------------------------

def _job_name(inh: str, mut: str, rep: int, chunk: int | None = None) -> str:
    base = f"{inh}_{mut}_r{rep}"
    return f"{base}_c{chunk}" if chunk is not None else base


def _amber_setup_line(amber_home: str) -> str:
    if amber_home:
        return f'source "{amber_home}/amber.sh"'
    return "# pmemd.cuda assumed to be on PATH"


def _write_run_scripts(
    replica_dir: Path,
    inh: str, mut: str, rep: int,
    sim: dict, slurm: dict, topology: str,
    mode: str,
    has_ligand: bool = True,
) -> None:
    use_hmr  = bool(sim.get("use_hmr", True))
    has_nmr  = bool(_get_nmr_list(sim, inh, has_ligand))

    hmr_block = (
        "log '00_prep: HMR'\ncpptraj -i HMR.ccptraj"
        if use_hmr else "# HMR disabled"
    )
    restrainer_block = (
        'cd "$DIR"\nlog "Generating NMR restraints"\npython restrainer.py\n'
        if has_nmr else ""
    )

    prod   = sim["production"]
    total  = int(prod["total_chunks"])
    cpj    = int(prod["chunks_per_job"])
    n_jobs = math.ceil(total / cpj)
    nvt_dir = replica_dir / "04_NVT"

    if mode == "cluster":
        module   = slurm["amber_module"]
        slurm_gpu = slurm["gpu"]
        account  = slurm["account"]

        header = fill(RUN_GPU_HEADER,
            WALLTIME=slurm_gpu["time"], NTASKS=slurm_gpu["ntasks"],
            GRES=slurm_gpu["gres"], PARTITION=slurm_gpu["partition"],
            JOBNAME=_job_name(inh, mut, rep),
            ACCOUNT=account, MODULE=module,
            REPLICA_DIR=str(replica_dir.resolve()),
        )
        body = fill(RUN_BODY,
            HMR_BLOCK=hmr_block,
            RESTRAINER_BLOCK=restrainer_block,
            MAXCAP=sim["min"]["max_cycles_cap"],
            CONVTHRESH=sim["min"]["convergence_threshold"],
            TOPOLOGY=topology,
            NVT_LAUNCH="sbatch run_NVT_1.cmd",
        )
        _write_exe(replica_dir / "run_gpu", header + body)

        for job_idx in range(1, n_jobs + 1):
            start    = (job_idx - 1) * cpj + 1
            end      = min(job_idx * cpj, total)
            next_job = f"\nsbatch run_NVT_{job_idx + 1}.cmd" if end < total else ""
            nvt_header = fill(NVT_CLUSTER_HEADER,
                WALLTIME=slurm_gpu["time"], NTASKS=slurm_gpu["ntasks"],
                GRES=slurm_gpu["gres"], PARTITION=slurm_gpu["partition"],
                JOBNAME=_job_name(inh, mut, rep, chunk=job_idx),
                ACCOUNT=account, MODULE=module,
            )
            nvt_body = fill(NVT_JOB_BODY,
                START=start, END=end, TOTAL=total,
                TOPOLOGY=topology, NEXTJOB=next_job,
            )
            _write_exe(nvt_dir / f"run_NVT_{job_idx}.cmd", nvt_header + nvt_body)

    else:  # local
        amber_setup = _amber_setup_line(sim.get("amber_home", ""))

        header = fill(RUN_LOCAL_HEADER, AMBER_SETUP=amber_setup)
        body   = fill(RUN_BODY,
            HMR_BLOCK=hmr_block,
            RESTRAINER_BLOCK=restrainer_block,
            MAXCAP=sim["min"]["max_cycles_cap"],
            CONVTHRESH=sim["min"]["convergence_threshold"],
            TOPOLOGY=topology,
            NVT_LAUNCH='bash "$DIR/04_NVT/run_NVT_1_local.sh"',
        )
        _write_exe(replica_dir / "run_local", header + body)

        for job_idx in range(1, n_jobs + 1):
            start    = (job_idx - 1) * cpj + 1
            end      = min(job_idx * cpj, total)
            next_job = (f'\nbash "$DIR/run_NVT_{job_idx + 1}_local.sh"'
                        if end < total else "")
            nvt_header = fill(NVT_LOCAL_HEADER,
                AMBER_SETUP=amber_setup,
                START=start, END=end, TOTAL=total,
            )
            nvt_body = fill(NVT_JOB_BODY,
                START=start, END=end, TOTAL=total,
                TOPOLOGY=topology, NEXTJOB=next_job,
            )
            _write_exe(nvt_dir / f"run_NVT_{job_idx}_local.sh", nvt_header + nvt_body)


def _write_exe(path: Path, content: str) -> None:
    path.write_text(content)
    path.chmod(0o755)


# ---------------------------------------------------------------------------
# Per-replica setup
# ---------------------------------------------------------------------------

def setup_replica(
    inh:             str,
    mut:             str,
    rep:             int,
    sim:             dict,
    slurm:           dict,
    protein_dir:     Path,
    ligands_lib_dir: Path,
    simulations_dir: Path,
    mode:            str = "cluster",
) -> None:
    replica_dir = simulations_dir / inh / mut / f"replica_{rep}"
    for d in _SUBDIRS:
        (replica_dir / d).mkdir(parents=True, exist_ok=True)

    # PDB -copy directly from protein/ (already clean from senda-complex)
    src_pdb = protein_dir / f"{mut}_{inh}_dimer.pdb"
    if not src_pdb.exists():
        raise FileNotFoundError(f"Protein PDB not found: {src_pdb}")
    shutil.copy(src_pdb, replica_dir / "00_prep" / src_pdb.name)

    # Ligand FF parameters (not needed for APO)
    has_ligand = (inh != "APO")
    if has_ligand:
        for ext in (".lib", ".frcmod"):
            src = ligands_lib_dir / inh / f"{inh}{ext}"
            if not src.exists():
                raise FileNotFoundError(
                    f"Ligand parameter file not found: {src}\n"
                    f"Run senda-param first to generate GAFF parameters."
                )
            shutil.copy(src, replica_dir / "00_prep" / src.name)

    # tleap input
    leap_cfg = sim.get("leap") or {}
    (replica_dir / "00_prep" / "leap_structure").write_text(
        _render_leap(inh, mut, has_ligand, leap_cfg)
    )

    # HMR cpptraj script
    if sim.get("use_hmr", True):
        (replica_dir / "00_prep" / "HMR.ccptraj").write_text(HMR_CCPTRAJ)

    # NMR restrainer (runs at job time, after tleap builds the topology)
    nmr_list = _get_nmr_list(sim, inh, has_ligand)
    if nmr_list:
        _write_exe(replica_dir / "restrainer.py", make_restrainer(nmr_list))

    # AMBER input files
    topology = _write_input_files(replica_dir, inh, sim, has_ligand=has_ligand)

    # Run scripts
    _write_run_scripts(replica_dir, inh, mut, rep, sim, slurm, topology, mode,
                       has_ligand=has_ligand)


# ---------------------------------------------------------------------------
# Batch setup and submit
# ---------------------------------------------------------------------------

def _master_script(replica_dir: Path, mode: str) -> Path:
    return replica_dir / ("run_gpu" if mode == "cluster" else "run_local")


def setup_all(
    inhibitors:      list[str],
    mutants:         list[str],
    n_replicas:      int,
    sim:             dict,
    slurm:           dict,
    protein_dir:     Path,
    ligands_lib_dir: Path,
    simulations_dir: Path,
    mode:            str  = "cluster",
    force:           bool = False,
) -> None:
    total = len(inhibitors) * len(mutants) * n_replicas
    print(f"\nMode: {mode}")
    print(f"Setting up {len(inhibitors)} inhibitor(s) x {len(mutants)} mutant(s) "
          f"x {n_replicas} replica(s) = {total} replica directories\n")

    use_hmr = sim.get("use_hmr", True)
    dt      = 0.004 if use_hmr else 0.002
    prod    = sim["production"]

    for inh in inhibitors:
        for mut in mutants:
            for rep in range(1, n_replicas + 1):
                label  = f"  {inh}/{mut}/replica_{rep}"
                r_dir  = simulations_dir / inh / mut / f"replica_{rep}"
                script = _master_script(r_dir, mode)
                if not force and script.exists():
                    print(f"{label}  (already exists -skipping)")
                    continue
                print(label)
                setup_replica(
                    inh, mut, rep, sim, slurm,
                    protein_dir, ligands_lib_dir, simulations_dir,
                    mode=mode,
                )

    print("\nSetup complete.")
    print(f"  HMR: {use_hmr}  (dt = {dt} ps)")
    print(f"  Production: {prod['total_chunks']} x {prod['ns_per_chunk']} ns "
          f"= {prod['total_chunks'] * prod['ns_per_chunk']} ns per replica")


def submit_all(
    inhibitors:      list[str],
    mutants:         list[str],
    n_replicas:      int,
    simulations_dir: Path,
    mode:            str = "cluster",
) -> None:
    for inh in inhibitors:
        for mut in mutants:
            for rep in range(1, n_replicas + 1):
                r_dir  = simulations_dir / inh / mut / f"replica_{rep}"
                script = _master_script(r_dir, mode)
                label  = f"  {inh}/{mut}/replica_{rep}"

                if not script.exists():
                    name = script.name
                    print(f"{label}  SKIP  ({name} not found -run setup first)")
                    continue

                if mode == "cluster":
                    result = subprocess.run(
                        ["sbatch", "run_gpu"],
                        capture_output=True, text=True,
                        cwd=str(r_dir),
                    )
                    if result.returncode == 0:
                        print(f"{label}  ->{result.stdout.strip()}")
                    else:
                        print(f"{label}  ERROR: {result.stderr.strip()}")
                else:
                    log_file = r_dir / "run.log"
                    proc = subprocess.Popen(
                        ["bash", str(script)],
                        stdout=log_file.open("w"),
                        stderr=subprocess.STDOUT,
                        start_new_session=True,
                    )
                    print(f"{label}  ->PID {proc.pid}  (log: {log_file})")

"""
senda.qmmm.string.prod

Stage 05_QMMM_restraint_free: optional unrestrained QM/MM production run.

Generates simulations/{inh}/{mut}/05_QMMM_restraint_free/ containing:
  in      -- AMBER QM/MM input (no restraints, restart from equil)
  prod.sh -- SLURM job script

Reads QM region metadata written by stage 05 (equil), so that stage must
run first.  If this stage output exists, the scan stage (06) will use
prod.rst7 as its starting structure instead of 0e.rst7.
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .equil import _load_stage_metadata
from ..common.templates import fill, sbatch_lines, STAGE05_IN, PROD_SLURM


def setup(
    inh:      str,
    mut:      str,
    inh_cfg:  dict,
    cfg:      dict,
    cwd:      Path,
    submit:   bool = False,
    after:    str | None = None,
) -> None:
    """
    Set up the restraint-free QM/MM production stage for one inhibitor/mutant pair.

    Reads metadata from stage 05 (equil).  No atom resolution is repeated here.
    """
    sim_base = cwd / "simulations" / inh / mut
    meta     = _load_stage_metadata(sim_base)
    qmmask   = meta["qmmask"]
    qmcharge = meta["qmcharge"]

    stage_dir = sim_base / "05_QMMM_restraint_free"
    stage_dir.mkdir(parents=True, exist_ok=True)

    string_cfg_top = (cfg.get("qmmm") or {}).get("string") or {}
    prod_cfg  = inh_cfg.get("prod")  or {}
    equil_cfg = inh_cfg.get("equil") or string_cfg_top.get("equil") or {}

    in_text = fill(
        STAGE05_IN,
        IREST     = 1,
        NTX       = 5,
        NMROPT    = "",
        DISANG    = "",
        TEMP      = prod_cfg.get("temp",     equil_cfg.get("temp",     300.0)),
        QMCUT     = inh_cfg.get("qmcut") or string_cfg_top.get("qmcut") or 12.0,
        GAMMA_LN  = prod_cfg.get("gamma_ln", 1.0),
        NSTLIM    = prod_cfg.get("nstlim",   100000),
        DT        = prod_cfg.get("dt",        0.001),
        NTPR      = prod_cfg.get("ntpr",      500),
        NTWX      = prod_cfg.get("ntwx",      500),
        NTWR      = prod_cfg.get("ntwr",      500),
        QMMASK    = qmmask,
        QMCHARGE  = qmcharge,
        QM_THEORY = meta["qm_theory"],
    )
    (stage_dir / "in").write_text(in_text)

    slurm_cfg  = cfg.get("slurm") or {}
    qmmm_cfg   = slurm_cfg.get("qmmm") or {}
    cpu_cfg    = slurm_cfg.get("cpu")  or {}

    rel_parm = f"../replica_1/00_prep/{meta['top_name']}"

    slurm_text = fill(
        PROD_SLURM,
        TIME          = cpu_cfg.get("time",   "3-00:00:00"),
        SCHEME        = meta["scheme"],
        NTASKS        = qmmm_cfg.get("ntasks", 8),
        EXTRA_SBATCH  = sbatch_lines(cpu_cfg, account=slurm_cfg.get("account")),
        ENV_SETUP     = slurm_cfg.get("env_setup", ""),
        PARM          = rel_parm,
    )
    script = stage_dir / "prod.sh"
    script.write_text(slurm_text)
    script.chmod(0o755)

    print(f"  Stage 05_QMMM_restraint_free written: {stage_dir}")

    if submit:
        dep = ["--dependency", f"afterok:{after}"] if after else []
        result = subprocess.run(
            ["sbatch"] + dep + [str(script)],
            capture_output=True, text=True, cwd=stage_dir,
        )
        print(f"  sbatch: {result.stdout.strip() or result.stderr.strip()}")

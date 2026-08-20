"""
senda.qmmm.string.string

Stage 07: adaptive string method QM/MM simulation.

Generates simulations/{inh}/{mut}/07_QMMM_string/ containing:
  in        -- AMBER string input (seed placeholder resolved by in.sh)
  in.sh     -- bash script: writes per-node {i}.in files and string.groupfile
  guess     -- string guess file (with AMBER header: N_nodes  N_cvs  0.0)
  CVs       -- AMBER CVs file (collective variables for sander ASM)
  string.sh -- SLURM job script
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from .equil import _load_stage_metadata
from ..common.cvs import write_string_guess, write_cvs_file
from ..common.templates import fill, sbatch_lines, DEFAULT_ENV_SETUP_STRING, STRING_IN, STRING_IN_SH, STRING_SLURM


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
    Set up the adaptive string method stage for one inhibitor/mutant pair.

    Reads metadata and cached guess from stage 05.
    Writes the AMBER input, groupfile generator, guess, CVs, and SLURM script.
    """
    import numpy as np

    sim_base  = cwd / "simulations" / inh / mut
    meta      = _load_stage_metadata(sim_base)
    guess     = np.load(str(sim_base / "_guess_interpolated.npy"))  # (n_nodes, n_cvs)

    n_nodes           = meta["n_nodes"]
    cv_indices_per_cv = meta["cv_indices_per_cv"]
    cv_specs          = inh_cfg.get("collective_variables") or []

    if guess.shape[0] != n_nodes:
        raise ValueError(
            f"Cached guess has {guess.shape[0]} rows but n_nodes={n_nodes}. "
            "Re-run stage 05 (equil) to rebuild the cache."
        )

    stage_dir = sim_base / "07_QMMM_string"
    stage_dir.mkdir(parents=True, exist_ok=True)

    string_cfg_top = (cfg.get("qmmm") or {}).get("string") or {}
    string_cfg = inh_cfg.get("string") or string_cfg_top.get("string") or {}
    equil_cfg  = inh_cfg.get("equil")  or string_cfg_top.get("equil")  or {}
    slurm_cfg  = cfg.get("slurm") or {}
    qmmm_cfg   = slurm_cfg.get("qmmm") or {}
    cpu_cfg    = slurm_cfg.get("cpu")  or {}

    # AMBER string input (seed resolved per-node by in.sh)
    seed = int(string_cfg.get("seed", 1234))
    in_text = fill(
        STRING_IN,
        TEMP             = string_cfg.get("temp",             equil_cfg.get("temp", 300.0)),
        QMCUT            = inh_cfg.get("qmcut") or string_cfg_top.get("qmcut") or 12.0,
        GAMMA_LN         = string_cfg.get("gamma_ln",         equil_cfg.get("gamma_ln", 5.0)),
        NSTLIM           = string_cfg.get("nstlim",           50000),
        DT               = string_cfg.get("dt",               0.001),
        NTPR             = string_cfg.get("ntpr",             100),
        NTWX             = string_cfg.get("ntwx",             100),
        NTWR             = string_cfg.get("ntwr",             100),
        QMMASK           = meta["qmmask"],
        QMCHARGE         = meta["qmcharge"],
        QM_THEORY        = meta["qm_theory"],
        PREP_STEPS       = string_cfg.get("prep_steps",       500),
        Z_BIAS           = str(string_cfg.get("z_bias", "false")).lower(),
        FORCE_CONSTANT_D = string_cfg.get("force_constant_d", 100.0),
    )
    (stage_dir / "in").write_text(in_text)

    # in.sh: writes per-node .in files and string.groupfile
    rel_parm_h10 = f"../{meta['h10_name']}"
    in_sh_text = fill(
        STRING_IN_SH,
        PARM_H10 = rel_parm_h10,
        SEED     = seed,
        N_NODES  = n_nodes,
    )
    in_sh = stage_dir / "in.sh"
    in_sh.write_text(in_sh_text)
    in_sh.chmod(0o755)

    # Guess file (with AMBER string header)
    write_string_guess(guess, stage_dir / "guess")

    # CVs file
    write_cvs_file(cv_specs, cv_indices_per_cv, stage_dir / "CVs")

    # SLURM script
    # ntasks for string = n_nodes * ntasks_per_node (2 MPI tasks per node by default)
    ntasks_per_node = int(qmmm_cfg.get("ntasks", 2))
    ntasks_string   = n_nodes * ntasks_per_node

    time_string = cpu_cfg.get("time_string") or cpu_cfg.get("time", "7-00:00:00")

    slurm_text = fill(
        STRING_SLURM,
        TIME          = time_string,
        SCHEME        = meta["scheme"],
        NTASKS_STRING = ntasks_string,
        EXTRA_SBATCH  = sbatch_lines(cpu_cfg, account=slurm_cfg.get("account")),
        ENV_SETUP     = slurm_cfg.get("env_setup", DEFAULT_ENV_SETUP_STRING),
        N_NODES       = n_nodes,
    )
    script = stage_dir / "string.sh"
    script.write_text(slurm_text)
    script.chmod(0o755)

    print(f"  Stage 07 written: {stage_dir}")

    if submit:
        dep = ["--dependency", f"afterok:{after}"] if after else []
        result = subprocess.run(
            ["sbatch"] + dep + [str(script)],
            capture_output=True, text=True, cwd=stage_dir,
        )
        print(f"  sbatch: {result.stdout.strip() or result.stderr.strip()}")

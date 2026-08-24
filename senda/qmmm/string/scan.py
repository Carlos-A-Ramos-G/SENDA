"""
senda.qmmm.string.scan

Stage 06: restrained QM/MM scan along the reaction path.

Generates simulations/{inh}/{mut}/06_QMMM_scan/ containing:
  in_template  -- AMBER input with __NODE__ placeholder
  restr0       -- extra_restraints as AMBER &rst blocks (if any)
  restr{i}     -- per-node CV restraints (from guess) + restr0 appended in job
  scan.sh      -- SLURM job script (sequential over all nodes)
  center.sh    -- cpptraj centering for one node; scan.sh calls it on
                  node 0 before the loop starts, then after every node's
                  sander.MPI run -- each node starts from the previous
                  node's centered structure, not the raw one
"""
from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np

from .equil import _load_stage_metadata
from ..common.cvs import build_rst_block
from ..common.templates import fill, sbatch_lines, DEFAULT_ENV_SETUP, SCAN_IN_TEMPLATE, SCAN_SLURM, CENTER_SH


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
    Set up the restrained scan stage for one inhibitor/mutant pair.

    Reads metadata and cached guess from stage 05.
    Writes all restraint files, the AMBER input template, and the SLURM script.
    """
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

    stage_dir = sim_base / "06_QMMM_scan"
    stage_dir.mkdir(parents=True, exist_ok=True)

    string_cfg_top = (cfg.get("qmmm") or {}).get("string") or {}
    scan_cfg = inh_cfg.get("scan") or string_cfg_top.get("scan") or {}
    force_constant = float(scan_cfg.get("force_constant", 100.0))

    # AMBER input template (NODE filled by scan.sh via sed)
    equil_cfg = inh_cfg.get("equil") or string_cfg_top.get("equil") or {}
    in_text = fill(
        SCAN_IN_TEMPLATE,
        TEMP       = scan_cfg.get("temp",     equil_cfg.get("temp",     300.0)),
        QMCUT      = inh_cfg.get("qmcut") or string_cfg_top.get("qmcut") or 12.0,
        GAMMA_LN   = scan_cfg.get("gamma_ln", equil_cfg.get("gamma_ln", 5.0)),
        NSTLIM     = scan_cfg.get("nstlim",   5000),
        DT         = scan_cfg.get("dt",       0.001),
        NTPR       = scan_cfg.get("ntpr",     50),
        NTWX       = scan_cfg.get("ntwx",     100),
        NTWR       = scan_cfg.get("ntwr",     100),
        QMMASK     = meta["qmmask"],
        QMCHARGE   = meta["qmcharge"],
        QM_THEORY  = meta["qm_theory"],
    )
    (stage_dir / "in_template").write_text(in_text)

    # Per-node restraint files
    for node_i in range(1, n_nodes + 1):
        target_row = guess[node_i - 1]  # CV target values for this node
        blocks = []
        for cv_idx, (cv, indices, target) in enumerate(
            zip(cv_specs, cv_indices_per_cv, target_row)
        ):
            cv_type = cv.get("type", "distance").lower()
            blocks.append(build_rst_block(indices, target, cv_type, force_constant))
        (stage_dir / f"restr{node_i}").write_text("".join(blocks))

    # Extra restraints file (appended to each restr{i} by the scan job)
    extra = inh_cfg.get("extra_restraints") or []
    restr0_text = _build_extra_restr(extra)
    (stage_dir / "restr0").write_text(restr0_text)

    # SLURM script
    slurm_cfg  = cfg.get("slurm") or {}
    qmmm_cfg   = slurm_cfg.get("qmmm") or {}
    cpu_cfg    = slurm_cfg.get("cpu")  or {}

    rel_parm = f"../replica_1/00_prep/{meta['top_name']}"

    slurm_text = fill(
        SCAN_SLURM,
        TIME          = cpu_cfg.get("time",   "2-00:00:00"),
        SCHEME        = meta["scheme"],
        NTASKS        = qmmm_cfg.get("ntasks", 8),
        EXTRA_SBATCH  = sbatch_lines(cpu_cfg, account=slurm_cfg.get("account")),
        ENV_SETUP     = slurm_cfg.get("env_setup", DEFAULT_ENV_SETUP),
        N_NODES       = n_nodes,
        PARM          = rel_parm,
    )
    script = stage_dir / "scan.sh"
    script.write_text(slurm_text)
    script.chmod(0o755)

    # Center.sh (cpptraj centering, called once per node right after that
    # node's sander.MPI run finishes)
    center_text = fill(
        CENTER_SH,
        PARM             = rel_parm,
        PROTEIN_LAST_RES = meta["n_protein_res"],
    )
    csh = stage_dir / "center.sh"
    csh.write_text(center_text)
    csh.chmod(0o755)

    print(f"  Stage 06 written: {stage_dir}")

    if submit:
        dep = ["--dependency", f"afterok:{after}"] if after else []
        result = subprocess.run(
            ["sbatch"] + dep + [str(script)],
            capture_output=True, text=True, cwd=stage_dir,
        )
        print(f"  sbatch: {result.stdout.strip() or result.stderr.strip()}")


# ---------------------------------------------------------------------------
# Restraint file helpers
# ---------------------------------------------------------------------------

def _build_extra_restr(extra_restraints: list) -> str:
    if not extra_restraints:
        return ""
    lines = []
    for r in extra_restraints:
        atoms = " ".join(str(a) for a in r["atoms"])
        lines.append("&rst")
        lines.append(f" iat={atoms},")
        lines.append(
            f" r1={r['r1']:.4f}, r2={r['r2']:.4f}, "
            f"r3={r['r3']:.4f}, r4={r['r4']:.4f},"
        )
        lines.append(f" rk2={r['rk2']:.2f}, rk3={r['rk3']:.2f},")
        lines.append("/")
    return "\n".join(lines) + "\n"

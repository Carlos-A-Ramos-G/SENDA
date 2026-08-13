"""
senda.qmmm.string.equil

Stage 05: QM/MM equilibration setup.

Generates simulations/{inh}/{mut}/05_QMMM_equilibration/ containing:
  in                -- AMBER QM/MM input
  restr             -- extra_restraints, plus soft CV restraints targeting
                       guess node 1 if equil.restrain_cvs is set (if any)
  equilibration.sh  -- SLURM job script

Also writes (at the mutant level):
  structure_H10.parm7  -- built from replica_1/00_prep/structure.parm7 (never
                          structure_HMR.parm7) with the reactive hydrogens'
                          mass set to 10 amu, for the string method's
                          mass-weighted path CV
"""
from __future__ import annotations

import subprocess
from pathlib import Path

from ..common.atoms import (
    resolve_cv_atoms_flat,
    resolve_cv_atoms_per_cv,
    resolve_qm_region,
    find_h10_atoms,
    count_protein_residues,
)
from ..common.h10 import generate_h10_topology
from ..common.cvs import (
    load_guess, interpolate_guess, write_scan_guess, write_string_guess,
    write_cvs_file, build_rst_block,
)
from ..common.templates import fill, STAGE05_IN, EQUIL_SLURM


def setup(
    inh:      str,
    mut:      str,
    inh_cfg:  dict,
    cfg:      dict,
    cwd:      Path,
    submit:   bool = False,
) -> None:
    """
    Set up the QM/MM equilibration stage for one inhibitor/mutant pair.

    Resolves atoms, builds the H10 topology, writes all input files,
    and optionally submits the SLURM job.
    """
    from senda.analysis.distances import (
        _find_pdb,
        _build_chain_map, _build_substrate_chain_map,
    )
    import pytraj as pt

    sim_base    = cwd / "simulations" / inh / mut
    rep1_dir    = sim_base / "replica_1"
    mc_cfg      = cfg.get("michaelis_complex") or {}
    protein_dir = cwd / mc_cfg.get("output_dir", "protein")
    chains      = mc_cfg.get("chains") or ["A"]
    chain       = chains[0]  # QM/MM uses chain A representative

    # Always build from the non-HMR topology, never structure_HMR.parm7:
    # the string method's mass-weighted path CV relies on H10 setting the
    # reactive hydrogens' mass to 10 amu relative to the true 1.008 amu
    # baseline. Starting from HMR's already-redistributed hydrogen masses
    # would miscalibrate that weighting for every other H in the system.
    top_path = rep1_dir / "00_prep" / "structure.parm7"
    if not top_path.exists():
        raise FileNotFoundError(f"Topology not found: {top_path}")
    pdb_path    = _find_pdb(rep1_dir, protein_dir, inh, mut)
    top         = pt.load_topology(str(top_path))
    top_atoms   = list(top.atoms)
    top_residues = list(top.residues)
    chain_map   = _build_chain_map(pdb_path, top)

    # Representative rst7 from senda-analyse output
    rst7_path = sim_base / f"{inh}_{mut}_chain{chain}_representative.rst7"
    if not rst7_path.exists():
        raise FileNotFoundError(
            f"Representative structure not found: {rst7_path}\n"
            "Run senda-analyse first."
        )

    rst7_traj  = pt.load(str(rst7_path), top=str(top_path))
    rst7_coords = rst7_traj.xyz[0]  # (n_atoms, 3)

    cv_specs = inh_cfg.get("collective_variables") or []
    if not cv_specs:
        raise ValueError(f"No collective_variables defined for {inh} in qmmm.string")

    # Collect substrate_sequence specs for proximity assignment
    seen_seq:       set  = set()
    sequence_specs: list = []
    substrate_A_set: set = set()

    def _collect(atom_specs):
        for spec in atom_specs:
            if "nearest_water_to" in spec:
                _collect([spec["nearest_water_to"]])
            elif "substrate_sequence" in spec:
                substrate_A_set.add(spec["substrate_sequence"])
            elif "sequence" in spec:
                key = (spec["sequence"], spec["name"])
                if key not in seen_seq:
                    seen_seq.add(key)
                    sequence_specs.append({"pdb_resnum": spec["sequence"], "name": spec["name"]})

    for cv in cv_specs:
        _collect(cv["atoms"])

    substrate_chain_map: dict = {}
    if substrate_A_set:
        sub_name            = top_residues[next(iter(substrate_A_set)) - 1].name
        all_substrate_amber = {r.index + 1 for r in top_residues if r.name == sub_name}
        substrate_chain_map = _build_substrate_chain_map(
            substrate_A_set, all_substrate_amber, sequence_specs,
            chain_map, chains, top_atoms, top_path, sim_base,
            coords=rst7_coords,
        )

    # Resolve CV atom indices (chain A)
    cv_indices_per_cv = resolve_cv_atoms_per_cv(
        cv_specs, chain_map, top_atoms, top_residues,
        substrate_chain_map, chain, rst7_coords,
    )
    cv_indices_flat = [i for per_cv in cv_indices_per_cv for i in per_cv]

    # Print resolved CV atoms for verification
    print(f"\n  CV atoms resolved for {inh}/{mut} chain {chain}:")
    for cv, indices in zip(cv_specs, cv_indices_per_cv):
        label = cv.get("type", "distance")
        names = [f"@{i}({top_atoms[i-1].name})" for i in indices]
        print(f"    {label}: {' -- '.join(names)}")

    # QM region and charge
    qmmask, qmcharge = resolve_qm_region(
        inh_cfg, cv_indices_flat, top_path,
        chain_map=chain_map, top_atoms=top_atoms, top_residues=top_residues,
        substrate_chain_map=substrate_chain_map, chain=chain, rst7_coords=rst7_coords,
    )

    # H10 topology
    h10_atoms = find_h10_atoms(cv_indices_flat, top_atoms, top_residues)
    h10_path  = sim_base / "structure_H10.parm7"
    generate_h10_topology(top_path, h10_atoms, h10_path)

    # Guess interpolation (write here so scan can reuse)
    string_cfg_top = (cfg.get("qmmm") or {}).get("string") or {}
    string_cfg = inh_cfg.get("string") or string_cfg_top.get("string") or {}
    scan_cfg   = inh_cfg.get("scan")   or string_cfg_top.get("scan")   or {}
    n_nodes    = int(scan_cfg.get("n_nodes", string_cfg.get("n_nodes", 64)))
    guess_path = Path(inh_cfg["guess"]) if not Path(inh_cfg.get("guess", "")).is_absolute() \
        else Path(inh_cfg["guess"])
    if not guess_path.is_absolute():
        guess_path = cwd / inh_cfg["guess"]
    guess_data = load_guess(guess_path)
    if guess_data.shape[0] != n_nodes:
        print(f"  Interpolating guess: {guess_data.shape[0]} -> {n_nodes} nodes")
        guess_data = interpolate_guess(guess_data, n_nodes)
    # Cache interpolated guess for scan/string stages
    guess_cache = sim_base / "_guess_interpolated.npy"
    import numpy as np
    np.save(str(guess_cache), guess_data)

    # Write CVs file (used by string stage but resolved here)
    n_protein_res = count_protein_residues(top_residues)

    # Create stage directory
    stage_dir = sim_base / "05_QMMM_equilibration"
    stage_dir.mkdir(parents=True, exist_ok=True)

    equil_cfg = inh_cfg.get("equil") or string_cfg_top.get("equil") or {}

    # Restraints (restr file): extra_restraints (absolute atom indices), plus
    # optional soft CV restraints targeting guess node 1 (the reactant-state
    # geometry) when equil.restrain_cvs is set.
    extra = inh_cfg.get("extra_restraints") or []
    restr_blocks = [_build_restr_file(extra)] if extra else []

    if equil_cfg.get("restrain_cvs"):
        force_constant = float(equil_cfg.get("force_constant", 20.0))
        cv_blocks = [
            build_rst_block(indices, float(target), cv.get("type", "distance").lower(), force_constant)
            for cv, indices, target in zip(cv_specs, cv_indices_per_cv, guess_data[0])
        ]
        restr_blocks.append("".join(cv_blocks))

    restr_text = "".join(restr_blocks)
    (stage_dir / "restr").write_text(restr_text)
    has_restraints = bool(restr_text.strip())

    # AMBER input
    in_text = fill(
        STAGE05_IN,
        IREST      = 0,
        NTX        = 1,
        NMROPT     = "\n  nmropt   = 1," if has_restraints else "",
        DISANG     = "&wt type = 'END'/\nDISANG=restr\n/\n" if has_restraints else "",
        TEMP       = equil_cfg.get("temp",     300.0),
        QMCUT      = inh_cfg.get("qmcut",      12.0),
        GAMMA_LN   = equil_cfg.get("gamma_ln", 5.0),
        NSTLIM     = equil_cfg.get("nstlim",   20000),
        DT         = equil_cfg.get("dt",        0.001),
        NTPR       = equil_cfg.get("ntpr",      50),
        NTWX       = equil_cfg.get("ntwx",      100),
        NTWR       = equil_cfg.get("ntwr",      100),
        QMMASK     = qmmask,
        QMCHARGE   = qmcharge,
        QM_THEORY  = inh_cfg.get("qm_theory",  "DFTB3"),
    )
    (stage_dir / "in").write_text(in_text)

    # SLURM script
    slurm_cfg  = cfg.get("slurm") or {}
    qmmm_slurm = slurm_cfg.get("qmmm") or {}
    cpu_cfg    = slurm_cfg.get("cpu")  or {}
    rel_rst7   = f"../{rst7_path.name}"
    rel_parm   = f"../replica_1/00_prep/{top_path.name}"

    senda_env = slurm_cfg.get("senda_env") or ""
    env_line  = senda_env if senda_env else "# (no senda_env set)"

    slurm_text = fill(
        EQUIL_SLURM,
        TIME         = qmmm_slurm.get("time",      "1-00:00:00"),
        SCHEME       = f"{inh}_{mut}",
        NTASKS       = qmmm_slurm.get("ntasks",    8),
        ACCOUNT      = slurm_cfg.get("account",    ""),
        PARTITION    = cpu_cfg.get("partition",     "cpu"),
        SENDA_ENV    = env_line,
        AMBER_MODULE = slurm_cfg.get("amber_module", "amber"),
        REP_RST7     = rel_rst7,
        PARM         = rel_parm,
    )
    script = stage_dir / "equilibration.sh"
    script.write_text(slurm_text)
    script.chmod(0o755)

    print(f"  Stage 05 written: {stage_dir}")

    if submit:
        result = subprocess.run(
            ["sbatch", str(script)], capture_output=True, text=True, cwd=stage_dir
        )
        print(f"  sbatch: {result.stdout.strip() or result.stderr.strip()}")

    # Store resolved metadata for downstream stages
    _write_stage_metadata(sim_base, {
        "qmmask":            qmmask,
        "qmcharge":          qmcharge,
        "qm_theory":         inh_cfg.get("qm_theory", "DFTB3"),
        "qmcut":             inh_cfg.get("qmcut", 12.0),
        "n_nodes":           n_nodes,
        "n_protein_res":     n_protein_res,
        "top_name":          top_path.name,
        "h10_name":          h10_path.name,
        "scheme":            f"{inh}_{mut}",
        "cv_indices_per_cv": cv_indices_per_cv,
    })


def _build_restr_file(extra_restraints: list) -> str:
    if not extra_restraints:
        return ""
    lines = []
    for r in extra_restraints:
        atoms = " ".join(str(a) for a in r["atoms"])
        lines.append("&rst")
        lines.append(f"iat={atoms}")
        lines.append(
            f"r1={r['r1']:.2f}, r2={r['r2']:.2f}, "
            f"r3={r['r3']:.2f}, r4={r['r4']:.2f}"
        )
        lines.append(f"rk2={r['rk2']}, rk3={r['rk3']}")
        lines.append("/")
    return "\n".join(lines) + "\n"


def _write_stage_metadata(sim_base: Path, meta: dict) -> None:
    import json
    path = sim_base / "_qmmm_string_meta.json"
    path.write_text(json.dumps(meta, indent=2))


def _load_stage_metadata(sim_base: Path) -> dict:
    import json
    path = sim_base / "_qmmm_string_meta.json"
    if not path.exists():
        raise FileNotFoundError(
            f"Stage metadata not found at {path}\n"
            "Run 'senda-qmmm string equil' first."
        )
    return json.loads(path.read_text())

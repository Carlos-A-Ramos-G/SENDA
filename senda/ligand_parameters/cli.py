"""
ligand_parameters.cli
Console-script entry points registered in pyproject.toml.

These commands handle individual ligands and implement the RESP workflow:

  ligand_parameters-run       [--ligand LIG] [--config config.yaml]
      Phase 1: write geometry-opt .com + resp.in + resp.qin for one ligand.

  ligand_parameters-hf-input  <opt.log> [--config config.yaml] [-c CHARGE]
      Phase 2b: parse Gaussian opt log -> write HF/6-31G(d) single-point .com.

  ligand_parameters-amber     <hf.log> [--config config.yaml] [--ligand LIG] [-c CHARGE]
      Phase 3: espgen -> resp -> antechamber -> remap names -> mol2 -> parmchk2 -> tleap.

  ligand_parameters-slurm     --ligand LIG [--config config.yaml]
      Generate a full-pipeline SLURM script for one ligand.

  ligand_parameters-check     [--config config.yaml]
      Verify that Gaussian and AmberTools are available.
"""

import sys


def _load_config(path):
    """Load config.yaml, or exit with a clear message if not found."""
    import yaml
    from pathlib import Path
    p = Path(path)
    if not p.exists():
        sys.exit(
            f"Error: config file '{path}' not found.\n"
            f"\n"
            f"  All ligand_parameters-* commands must be run from the directory that\n"
            f"  contains your config.yaml, or you must supply an explicit path:\n"
            f"\n"
            f"    ligand_parameters-run      --ligand LIG --config /path/to/config.yaml\n"
            f"    ligand_parameters-hf-input <opt.log>    --config /path/to/config.yaml\n"
            f"    ligand_parameters-amber    <hf.log>      --config /path/to/config.yaml\n"
            f"    ligand_parameters-slurm    --ligand LIG --config /path/to/config.yaml\n"
            f"\n"
            f"  Copy config.yaml from the ligand_parameters repository into\n"
            f"  your working directory and fill in the values for your system."
        )
    with open(p) as f:
        return yaml.safe_load(f)


# ---------------------------------------------------------------------------
# ligand_parameters-run   (Phase 1)
# ---------------------------------------------------------------------------

def run():
    """Write Gaussian geometry-opt input and RESP files for one ligand (Phase 1)."""
    import argparse
    from pathlib import Path
    from senda.ligand_parameters.gaussian import (write_com, NPROC_DEFAULT, MEM_DEFAULT,
                                       ROUTE_DEFAULT)
    from senda.ligand_parameters.resp import write_resp_in, write_resp_qin
    from senda.ligand_parameters._utils import get_resname, resolve_charge

    parser = argparse.ArgumentParser(
        description="Phase 1: write Gaussian opt input + RESP files for one ligand.",
    )
    parser.add_argument('--ligand', required=True,
                        help='Ligand ID (PDB filename stem, e.g. "LIG")')
    parser.add_argument('--config', default='config.yaml',
                        help='Config file (default: config.yaml)')
    parser.add_argument('-c', '--charge', type=int, default=None,
                        help='Override net charge (default: auto-detect via RDKit)')
    parser.add_argument('-m', '--mult', type=int, default=None,
                        help='Override spin multiplicity (default: from config)')
    args = parser.parse_args()

    cfg  = _load_config(args.config)
    base = Path(cfg.get('base_dir', '.')).resolve()
    lcfg = cfg.get('ligand', {}) or {}
    gcfg = cfg.get('gaussian_opt', {}) or cfg.get('gaussian', {}) or {}

    ligands_dir = base / cfg.get('ligands_dir', 'ligands')
    output_dir  = base / cfg.get('output_dir',  'output')
    ligand_pdb  = ligands_dir / f"{args.ligand}.pdb"

    if not ligand_pdb.exists():
        sys.exit(f"Error: ligand PDB not found: {ligand_pdb}")

    resname = get_resname(ligand_pdb)
    charge  = args.charge if args.charge is not None else resolve_charge(
        ligand_pdb, lcfg)
    mult    = args.mult   if args.mult   is not None else lcfg.get('multiplicity', 1)

    workdir = output_dir / args.ligand
    workdir.mkdir(parents=True, exist_ok=True)

    lig     = args.ligand
    opt_com = str(workdir / f"{lig}_opt.com")

    print("=" * 60)
    print(f"Ligand   : {lig}  (resname '{resname}',  charge {charge:+d})")
    print("=" * 60)

    print()
    print("=" * 60)
    print("Step 1 -- Geometry-optimisation Gaussian input")
    print("=" * 60)
    write_com(
        str(ligand_pdb), opt_com,
        charge=charge, mult=mult,
        nproc=gcfg.get('nproc', NPROC_DEFAULT),
        mem=gcfg.get('mem',   MEM_DEFAULT),
        route=gcfg.get('route', ROUTE_DEFAULT),
    )

    print()
    print("=" * 60)
    print("Steps 2-3 -- RESP input files")
    print("=" * 60)
    write_resp_in(str(ligand_pdb), charge, resname, str(workdir / 'resp.in'))
    write_resp_qin(str(ligand_pdb), str(workdir / 'resp.qin'))

    print()
    print(f"All Phase 1 outputs written to: {workdir}/")
    print()
    print("Next steps:")
    print(f"  1. Submit {workdir}/{lig}_opt.com to HPC")
    print(f"  2. Copy {lig}_opt.log into {workdir}/, then:")
    print(f"       ligand_parameters-hf-input {workdir}/{lig}_opt.log --config {args.config}")
    print(f"  3. Submit {workdir}/{lig}_hf.com to HPC")
    print(f"  4. Copy {lig}_hf.log into {workdir}/, then:")
    print(f"       ligand_parameters-amber {workdir}/{lig}_hf.log --config {args.config}")
    print()
    print(f"  Or generate a single end-to-end SLURM script:")
    print(f"       ligand_parameters-slurm --ligand {lig} --config {args.config}")


# ---------------------------------------------------------------------------
# ligand_parameters-hf-input   (Phase 2b)
# ---------------------------------------------------------------------------

def hf_input():
    """Extract optimised geometry from Gaussian log and write HF/ESP input."""
    import argparse
    from pathlib import Path
    from senda.ligand_parameters.gaussian import (parse_opt_log, write_hf_com,
                                       NPROC_DEFAULT, MEM_DEFAULT,
                                       HF_ROUTE_DEFAULT)

    parser = argparse.ArgumentParser(
        description="Gaussian opt log -> HF/6-31G(d) single-point .com",
    )
    parser.add_argument('log',            help='Gaussian optimisation log (e.g. LIG_opt.log)')
    parser.add_argument('com', nargs='?', help='Output .com  (default: <base>_hf.com)')
    parser.add_argument('--config',       default='config.yaml')
    parser.add_argument('-c', '--charge', type=int, default=None)
    parser.add_argument('-m', '--mult',   type=int, default=None)
    parser.add_argument('-n', '--nproc',  type=int, default=None)
    parser.add_argument('--mem',          default=None)
    args = parser.parse_args()

    cfg  = _load_config(args.config)
    lcfg = cfg.get('ligand', {}) or {}
    gcfg = cfg.get('gaussian_hf', {}) or cfg.get('gaussian', {}) or {}

    charge = args.charge if args.charge is not None else lcfg.get(
        'net_charge', 0 if str(lcfg.get('net_charge', 'auto')).lower() != 'auto' else 0)
    mult   = args.mult   if args.mult   is not None else lcfg.get('multiplicity', 1)
    nproc  = args.nproc  if args.nproc  is not None else gcfg.get('nproc', NPROC_DEFAULT)
    mem    = args.mem    if args.mem    is not None else gcfg.get('mem', MEM_DEFAULT)

    log_path = Path(args.log)
    workdir  = log_path.parent

    # Derive base name: strip _opt/_hf suffixes to recover the ligand stem
    stem = log_path.stem
    for tag in ('_opt', '_hf'):
        stem = stem.replace(tag, '')
    base     = stem
    com_path = args.com or str(workdir / f"{base}_hf.com")

    # If charge was 'auto', try to read from the config's resolved charge path
    if str(lcfg.get('net_charge', 'auto')).lower() == 'auto' and args.charge is None:
        cfg_base = Path(cfg.get('base_dir', '.')).resolve()
        ligands_dir = cfg_base / cfg.get('ligands_dir', 'ligands')
        ligand_pdb  = ligands_dir / f"{base}.pdb"
        if ligand_pdb.exists():
            from senda.ligand_parameters._utils import infer_net_charge
            try:
                charge = infer_net_charge(ligand_pdb)
                print(f"Auto-detected charge: {charge:+d} (from {ligand_pdb.name})")
            except RuntimeError as e:
                sys.exit(f"Error: {e}\nUse -c to specify the charge explicitly.")
        else:
            print(f"Warning: charge is 'auto' but {ligand_pdb} not found -- "
                  "using charge=0. Pass -c to override.")

    print(f"Parsing optimised geometry from {args.log} ...")
    atoms_xyz = parse_opt_log(args.log)
    print(f"  Found {len(atoms_xyz)} atoms in final Standard orientation block")

    print(f"Writing HF/6-31G(d) single-point input ...")
    write_hf_com(atoms_xyz, com_path, base,
                 charge=charge, mult=mult, nproc=nproc, mem=mem,
                 route=gcfg.get('route') or HF_ROUTE_DEFAULT)


# ---------------------------------------------------------------------------
# ligand_parameters-amber   (Phase 3)
# ---------------------------------------------------------------------------

def amber():
    """Run espgen -> resp -> antechamber -> remap -> mol2 -> parmchk2 -> tleap."""
    import argparse
    from pathlib import Path
    from senda.ligand_parameters.amber import run_amber_pipeline
    from senda.ligand_parameters._utils import get_resname, resolve_charge

    parser = argparse.ArgumentParser(
        description="Run AMBER parameterization pipeline from Gaussian HF log (Phase 3).",
    )
    parser.add_argument('log',             help='Gaussian HF/ESP log (e.g. LIG_hf.log)')
    parser.add_argument('--config',        default='config.yaml')
    parser.add_argument('--ligand',        default=None,
                        help='Ligand ID (default: derived from log filename stem)')
    parser.add_argument('-c', '--charge',  type=int, default=None,
                        help='Override net charge')
    parser.add_argument('--workdir',       default='.',
                        help='Working directory (default: same dir as log)')
    args = parser.parse_args()

    cfg  = _load_config(args.config)
    lcfg = cfg.get('ligand',     {}) or {}
    amb_cfg = cfg.get('amber',   {}) or {}
    ff_cfg  = cfg.get('forcefield', {}) or {}

    log_path = Path(args.log)

    # Resolve ligand ID from --ligand flag or log filename
    if args.ligand:
        lig = args.ligand
    else:
        stem = log_path.stem
        for tag in ('_hf', '_opt'):
            stem = stem.replace(tag, '')
        lig = stem

    cfg_base    = Path(cfg.get('base_dir', '.')).resolve()
    ligands_dir = cfg_base / cfg.get('ligands_dir', 'ligands')
    ligand_pdb  = ligands_dir / f"{lig}.pdb"

    if not ligand_pdb.exists():
        sys.exit(
            f"Error: ligand PDB not found: {ligand_pdb}\n"
            f"  Use --ligand to specify the ligand ID explicitly."
        )

    resname = get_resname(ligand_pdb)
    charge  = (args.charge if args.charge is not None
               else resolve_charge(ligand_pdb, lcfg))

    workdir = (args.workdir if args.workdir != '.'
               else str(log_path.resolve().parent))

    print(f"Ligand   : {lig}  (resname '{resname}',  charge {charge:+d})")
    print(f"Log      : {args.log}")
    print(f"Workdir  : {workdir}")
    print()

    run_amber_pipeline(
        hf_log              = args.log,
        resname             = resname,
        charge              = charge,
        ligand_pdb          = str(ligand_pdb),
        workdir             = workdir,
        atom_type           = lcfg.get('atom_type', 'gaff2'),
        forcefield          = ff_cfg.get('ligand', 'leaprc.gaff2'),
        atom_type_overrides = amb_cfg.get('atom_type_overrides') or {},
    )


# ---------------------------------------------------------------------------
# ligand_parameters-slurm
# ---------------------------------------------------------------------------

def slurm():
    """Generate a full-pipeline SLURM script for one ligand (RESP workflow)."""
    import argparse
    from pathlib import Path
    from senda.ligand_parameters.slurm import write_slurm

    parser = argparse.ArgumentParser(
        description="Generate an end-to-end SLURM script for the RESP parameterization "
                    "pipeline for one ligand.",
    )
    parser.add_argument('--ligand', required=True,
                        help='Ligand ID (PDB filename stem)')
    parser.add_argument('--config', default='config.yaml')
    args = parser.parse_args()

    config_path = Path(args.config).resolve()
    cfg  = _load_config(args.config)
    base = Path(cfg.get('base_dir', '.')).resolve()

    output_dir  = base / cfg.get('output_dir', 'output')
    scripts_dir = base / 'scripts'
    scripts_dir.mkdir(parents=True, exist_ok=True)

    workdir = output_dir / args.ligand
    workdir.mkdir(parents=True, exist_ok=True)

    output = scripts_dir / f"{args.ligand}_resp.sh"

    write_slurm(
        lig         = args.ligand,
        cfg         = cfg,
        output      = output,
        proj_root   = base,
        workdir     = workdir,
        config_path = config_path,
    )


# ---------------------------------------------------------------------------
# ligand_parameters-check
# ---------------------------------------------------------------------------

def check():
    """Check that all required tools and packages are available."""
    from senda.ligand_parameters.check import check_env
    check_env()

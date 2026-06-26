"""
senda.simulation.run
CLI entry point for senda-sim.

Usage:
  senda-sim setup   --config config.yaml                    # cluster mode (SLURM)
  senda-sim setup   --config config.yaml --mode local       # local GPU workstation
  senda-sim setup   --config config.yaml --submit           # generate + submit
  senda-sim setup   --config config.yaml --force            # overwrite existing
  senda-sim submit  --config config.yaml                    # submit existing scripts

The default mode is read from amber_simulator.execution_mode in config.yaml
(falls back to 'cluster' if absent).
"""

import logging
import sys
from pathlib import Path

log = logging.getLogger(__name__)


def _load_config(config_path: Path) -> dict:
    if not config_path.exists():
        sys.exit(f"Config not found: {config_path}")
    import yaml
    with open(config_path) as fh:
        return yaml.safe_load(fh) or {}


def main() -> None:
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", required=True, metavar="CONFIG",
        help="Path to the senda config.yaml")

    sub = parser.add_subparsers(dest="command", required=True)

    p_setup = sub.add_parser("setup", help="Generate all replica directories")
    p_setup.add_argument("--mode", choices=["cluster", "local"], default=None,
        help="cluster (SLURM) or local (single GPU workstation); "
             "default: amber_simulator.execution_mode from config, or 'cluster'")
    p_setup.add_argument("--submit", action="store_true",
        help="Launch jobs immediately after setup")
    p_setup.add_argument("--force", action="store_true",
        help="Overwrite existing replica directories")

    p_submit = sub.add_parser("submit", help="Launch existing run scripts")
    p_submit.add_argument("--mode", choices=["cluster", "local"], default=None,
        help="Must match the mode used during setup")

    args = parser.parse_args()
    cfg  = _load_config(Path(args.config))
    cwd  = Path.cwd()

    inhibitors = cfg.get("inhibitors") or []
    mutants    = cfg.get("mutants")    or []
    n_replicas = int(cfg.get("replicas", 1))

    if not inhibitors:
        sys.exit("Config error: 'inhibitors' list is empty or missing.")
    if not mutants:
        sys.exit("Config error: 'mutants' list is empty or missing.")

    sim   = cfg.get("amber_simulator") or {}
    slurm = cfg.get("slurm")           or {}

    # Resolve mode: CLI flag > config execution_mode > default cluster
    cfg_mode = sim.get("execution_mode", "cluster")
    mode     = args.mode or cfg_mode

    mc_cfg          = cfg.get("michaelis_complex") or {}
    protein_dir     = cwd / mc_cfg.get("output_dir", "protein")
    ligands_lib_dir = cwd / "ligands_libraries"
    simulations_dir = cwd / "simulations"

    print(f"Mode            : {mode}")
    print(f"Protein dir     : {protein_dir}")
    print(f"Ligand lib dir  : {ligands_lib_dir}")
    print(f"Simulations dir : {simulations_dir}")
    print(f"Inhibitors      : {inhibitors}")
    print(f"Mutants         : {mutants}")
    print(f"Replicas        : {n_replicas}")

    from .setup import setup_all, submit_all

    if args.command == "setup":
        setup_all(
            inhibitors, mutants, n_replicas,
            sim, slurm,
            protein_dir, ligands_lib_dir, simulations_dir,
            mode=mode, force=args.force,
        )
        if args.submit:
            submit_all(inhibitors, mutants, n_replicas, simulations_dir, mode=mode)

    elif args.command == "submit":
        submit_all(inhibitors, mutants, n_replicas, simulations_dir, mode=mode)

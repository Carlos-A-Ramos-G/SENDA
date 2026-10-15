"""
senda.simulation.openmm_engine.cli
Entry point for senda-omm, called by the generated run scripts on the compute node.

Usage:
  senda-omm run --stage min   --replica-dir <replica>
  senda-omm run --stage heat  --replica-dir <replica>
  senda-omm run --stage equil --replica-dir <replica>
  senda-omm run --stage prod  --chunk <i> --replica-dir <replica>

All parameters come from <replica>/openmm_params.json (written by
senda-sim setup), so the config YAML is not needed on the compute node.
Equivalent without the console script: python -m senda.simulation.openmm_engine
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

PARAMS_FILE = "openmm_params.json"


def main(argv: list[str] | None = None) -> None:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s")
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command", required=True)
    p_run = sub.add_parser("run", help="Run one MM stage")
    p_run.add_argument("--stage", required=True, choices=["min", "heat", "equil", "prod"])
    p_run.add_argument("--chunk", type=int, help="NVT production chunk (1-based, --stage prod)")
    p_run.add_argument("--replica-dir", type=Path, default=Path("."))
    args = parser.parse_args(argv)

    rd = args.replica_dir.resolve()
    params_path = rd / PARAMS_FILE
    if not params_path.exists():
        sys.exit(f"{params_path} not found -- run senda-sim setup with engine: openmm")
    params = json.loads(params_path.read_text())

    from . import stages
    if args.stage == "min":
        stages.minimize(params, rd)
    elif args.stage == "heat":
        stages.heat(params, rd)
    elif args.stage == "equil":
        stages.equilibrate(params, rd)
    else:
        if not args.chunk or args.chunk < 1:
            sys.exit("--stage prod requires --chunk >= 1")
        stages.produce(params, rd, args.chunk)


if __name__ == "__main__":
    main()

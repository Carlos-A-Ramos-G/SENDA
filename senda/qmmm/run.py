"""
senda.qmmm.run

Top-level entry point for `senda-qmmm`.

Usage:
  senda-qmmm string equil   --config config.yaml [options]
  senda-qmmm string scan    --config config.yaml [options]
  senda-qmmm string string  --config config.yaml [options]
"""
from __future__ import annotations

import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="senda-qmmm",
        description="QM/MM simulation workflow for senda",
    )
    sub = parser.add_subparsers(dest="method", required=True)
    sub.add_parser("string", help="Adaptive string method workflow",
                   add_help=False)

    # Parse only the first positional to dispatch; pass the rest through
    args, remaining = parser.parse_known_args(argv)

    if args.method == "string":
        from .string.run import main as string_main
        string_main(remaining)
    else:
        print(f"Unknown method: {args.method}", file=sys.stderr)
        sys.exit(1)

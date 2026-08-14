"""
senda.qmmm.run

Top-level entry point for `senda-qmmm`.

Usage:
  senda-qmmm --config config.yaml string equil   [options]
  senda-qmmm --config config.yaml string scan    [options]
  senda-qmmm --config config.yaml string string  [options]
"""
from __future__ import annotations

import argparse
import sys


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="senda-qmmm",
        description="QM/MM simulation workflow for senda",
    )
    parser.add_argument("--config", required=True, metavar="CONFIG",
                        help="Path to senda config.yaml")
    sub = parser.add_subparsers(dest="method", required=True)
    sub.add_parser("string", help="Adaptive string method workflow",
                   add_help=False)

    # Parse --config and the first positional to dispatch; pass the rest
    # through, re-injecting --config so the inner subparser sees it too.
    args, remaining = parser.parse_known_args(argv)

    if args.method == "string":
        from .string.run import main as string_main
        string_main(remaining + ["--config", args.config])
    else:
        print(f"Unknown method: {args.method}", file=sys.stderr)
        sys.exit(1)

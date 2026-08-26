"""
senda.integration.run

Entry point for `senda-integration` subcommands:
  string  -- PMF/free-energy integration of stage 07's (adaptive string
             method) sampling output: converts sampling/*_final.dat to
             WHAM/MBAR/vFEP input, runs chunked MBAR for a standard-error
             estimate, and plots the resulting PMF.

Usage:
  senda-integration --config config.yaml string [-i INH] [-m MUT]

Without -i/-m: processes every inhibitor under qmmm.string.inhibitors that's
also present in the top-level inhibitors: list (or all of them if that list
is empty), and every mutant from the per-inhibitor or top-level mutants:
list. -i/-m each explicitly select one inhibitor/mutant, bypassing those
filters entirely -- even for a pair not listed anywhere else in the config.

Requires the external ndfes/ndfes-PrintFES.py tools in $PATH, and
matplotlib (pip install -e ".[analysis]") for the PMF plot.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from senda.qmmm.string.run import _iter_pairs


def _load_config(config_path: Path) -> dict:
    from senda.config import resolve_slurm_profile
    with open(config_path) as fh:
        cfg = yaml.safe_load(fh)
    return resolve_slurm_profile(cfg)


def main_string(args):
    from .string import setup
    cfg = _load_config(args.config)
    cwd = args.config.parent
    for inh, mut, inh_cfg in _iter_pairs(cfg, args.inh, args.mut):
        print(f"\n[integration string] {inh}/{mut}")
        setup(inh, mut, inh_cfg, cfg, cwd)


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("-i", "--inh", default=None, metavar="INH",
                   help="Run only for this inhibitor -- bypasses the top-level "
                        "inhibitors: filter (default: every inhibitor under "
                        "qmmm.string.inhibitors that's also in the top-level "
                        "inhibitors: list, or all of them if that list is empty)")
    p.add_argument("-m", "--mut", default=None, metavar="MUT",
                   help="Run only for this mutant -- bypasses the per-inhibitor/"
                        "top-level mutants: fallback entirely, even for a mutant "
                        "not listed anywhere in the config")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="senda-integration",
        description="Free energy integration for QM/MM string method output",
    )
    parser.add_argument("--config", type=Path, required=True, metavar="CONFIG",
                        help="Path to senda config.yaml")
    sub = parser.add_subparsers(dest="stage", required=True)

    p_str = sub.add_parser("string", help="PMF integration of stage 07's sampling output")
    _add_common(p_str)

    args = parser.parse_args(argv)

    dispatch = {"string": main_string}
    try:
        dispatch[args.stage](args)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

"""
senda.qmmm.string.run

Entry point for `senda-qmmm` subcommands:
  equil   -- stage 05: QM/MM equilibration setup
  prod    -- stage 05_QMMM_restraint_free: optional unrestrained QM/MM production
  scan    -- stage 06: restrained scan setup
  string  -- stage 07: adaptive string method setup

Usage:
  senda-qmmm --config config.yaml equil  [-s] [-i INH] [-m MUT]
  senda-qmmm --config config.yaml prod   [-s] [-a JOBID] [-i INH] [-m MUT]
  senda-qmmm --config config.yaml scan   [-s] [-a JOBID] [-i INH] [-m MUT]
  senda-qmmm --config config.yaml string [-s] [-a JOBID] [-i INH] [-m MUT]

Without -i/-m: processes every inhibitor under qmmm.string.inhibitors that's
also present in the top-level inhibitors: list (or all of them if that list
is empty), and every mutant from the per-inhibitor or top-level mutants:
list. -i/-m each explicitly select one inhibitor/mutant, bypassing those
filters entirely -- even for a pair not listed anywhere else in the config.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml


def _load_config(config_path: Path) -> dict:
    with open(config_path) as fh:
        return yaml.safe_load(fh)


def _iter_pairs(cfg: dict, inh_filter: str | None, mut_filter: str | None):
    """Yield (inh, mut, inh_cfg) for all matching inhibitor/mutant pairs."""
    qmmm_cfg = cfg.get("qmmm") or {}
    string_cfg_top = qmmm_cfg.get("string") or {}
    inhibitors = string_cfg_top.get("inhibitors") or {}
    top_inhibitors = cfg.get("inhibitors") or []

    for inh, inh_cfg in inhibitors.items():
        # -i explicitly requests one inhibitor -- bypasses the top-level
        # inhibitors: gate entirely. Without -i, that gate applies (if set).
        if inh_filter:
            if inh != inh_filter:
                continue
        elif top_inhibitors and inh not in top_inhibitors:
            continue

        # -m explicitly requests one mutant -- bypasses the per-inhibitor/
        # top-level mutants: fallback chain entirely.
        if mut_filter:
            yield inh, mut_filter, inh_cfg
            continue

        mutants = inh_cfg.get("mutants") or cfg.get("mutants") or [None]
        for mut in mutants:
            if mut is None:
                mut = "WT"
            yield inh, mut, inh_cfg


def main_equil(args):
    from .equil import setup
    cfg     = _load_config(args.config)
    cwd     = args.config.parent
    for inh, mut, inh_cfg in _iter_pairs(cfg, args.inh, args.mut):
        print(f"\n[equil] {inh}/{mut}")
        setup(inh, mut, inh_cfg, cfg, cwd, submit=args.submit)


def main_prod(args):
    from .prod import setup
    cfg = _load_config(args.config)
    cwd = args.config.parent
    for inh, mut, inh_cfg in _iter_pairs(cfg, args.inh, args.mut):
        print(f"\n[prod] {inh}/{mut}")
        setup(inh, mut, inh_cfg, cfg, cwd, submit=args.submit, after=args.after)


def main_scan(args):
    from .scan import setup
    cfg = _load_config(args.config)
    cwd = args.config.parent
    for inh, mut, inh_cfg in _iter_pairs(cfg, args.inh, args.mut):
        print(f"\n[scan] {inh}/{mut}")
        setup(inh, mut, inh_cfg, cfg, cwd, submit=args.submit, after=args.after)


def main_string(args):
    from .string import setup
    cfg = _load_config(args.config)
    cwd = args.config.parent
    for inh, mut, inh_cfg in _iter_pairs(cfg, args.inh, args.mut):
        print(f"\n[string] {inh}/{mut}")
        setup(inh, mut, inh_cfg, cfg, cwd, submit=args.submit, after=args.after)


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
    p.add_argument("-s", "--submit", action="store_true",
                   help="Submit SLURM job after writing files")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="senda-qmmm",
        description="QM/MM adaptive string method setup subcommands",
    )
    parser.add_argument("--config", type=Path, required=True, metavar="CONFIG",
                        help="Path to senda config.yaml")
    sub = parser.add_subparsers(dest="stage", required=True)

    p_equil = sub.add_parser("equil",  help="Stage 05: QM/MM equilibration")
    _add_common(p_equil)

    p_prod  = sub.add_parser("prod",   help="Stage 05_QMMM_restraint_free: optional unrestrained production")
    _add_common(p_prod)
    p_prod.add_argument("-a", "--after", default=None, metavar="JOBID",
                        help="SLURM dependency: afterok:<JOBID>")

    p_scan  = sub.add_parser("scan",   help="Stage 06: restrained scan")
    _add_common(p_scan)
    p_scan.add_argument("-a", "--after", default=None, metavar="JOBID",
                        help="SLURM dependency: afterok:<JOBID>")

    p_str   = sub.add_parser("string", help="Stage 07: adaptive string method")
    _add_common(p_str)
    p_str.add_argument("-a", "--after", default=None, metavar="JOBID",
                       help="SLURM dependency: afterok:<JOBID>")

    args = parser.parse_args(argv)

    dispatch = {"equil": main_equil, "prod": main_prod, "scan": main_scan, "string": main_string}
    try:
        dispatch[args.stage](args)
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        sys.exit(1)

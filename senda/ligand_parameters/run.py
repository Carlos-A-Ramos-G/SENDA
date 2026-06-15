"""
senda.ligand_parameters.run
Batch ligand parameterization runner.

Discovers all .pdb files in ./ligands/ relative to the working directory and
parameterizes each with AMBER/GAFF. Charge method is controlled by config.yaml:

  charge_method: bcc   → AM1-BCC via antechamber (default, no HPC needed)
  charge_method: resp  → RESP charges from Gaussian QM (Phase 1 only: writes
                         Gaussian inputs; Gaussian itself must run on HPC)

Outputs go to ./ligands_libraries/<ligand>/.

Usage (from inside a project directory that contains a ligands/ folder):
    senda-param
    senda-param --ligands LER NIR
    senda-param --jobs 4
    senda-param --force
"""

import concurrent.futures
import logging
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

log = logging.getLogger(__name__)

_LIGANDS_DIR   = "ligands"
_OUTPUT_DIR    = "ligands_libraries"
_FORCEFIELD    = "leaprc.gaff2"
_CHARGE_METHOD = "bcc"
_MULTIPLICITY  = 1


# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

def _load_config(config_path: Path) -> tuple[dict, list[str] | None]:
    """
    Parse config.yaml and return (ligand_parameters section, inhibitors list).
    inhibitors is None when the key is absent from the config.
    """
    if not config_path.exists():
        sys.exit(f"Error: config file not found: {config_path}")
    import yaml
    with open(config_path) as f:
        raw = yaml.safe_load(f) or {}
    lp_cfg     = raw.get("ligand_parameters", {}) or {}
    inhibitors = raw.get("inhibitors") or None
    return lp_cfg, inhibitors


# ---------------------------------------------------------------------------
# Ligand discovery
# ---------------------------------------------------------------------------

def _discover_ligands(ligands_dir: Path) -> list[str]:
    pdbs = sorted(ligands_dir.glob("*.pdb"), key=lambda p: p.stem)
    if not pdbs:
        raise FileNotFoundError(f"No .pdb files found in {ligands_dir}")
    return [p.stem for p in pdbs]


# ---------------------------------------------------------------------------
# tleap input writer (BCC mode)
# ---------------------------------------------------------------------------

def _write_tleap_input(work: Path, lig: str, resname: str, forcefield: str) -> Path:
    mol2   = (work / f"{lig}.mol2").resolve()
    frcmod = (work / f"{lig}.frcmod").resolve()
    lib    = (work / f"{lig}.lib").resolve()
    tleap  = work / "tleap_params.in"
    tleap.write_text(
        f"source {forcefield}\n"
        f"{resname} = loadmol2 {mol2}\n"
        f"check {resname}\n"
        f"loadamberparams {frcmod}\n"
        f"saveoff {resname} {lib}\n"
        "quit\n"
    )
    return tleap


# ---------------------------------------------------------------------------
# RESP Phase 1: generate Gaussian inputs
# ---------------------------------------------------------------------------

def _run_resp_phase1(lig: str, cwd: Path, cfg: dict,
                     work: Path, ligand_pdb: Path,
                     charge: int, resname: str, mult: int) -> tuple[bool, str]:
    """
    Write the Gaussian geometry-opt .com, resp.in, and resp.qin for one ligand.
    Gaussian itself must be run on HPC — this step only generates the inputs.
    """
    from senda.ligand_parameters.gaussian import (
        write_com, NPROC_DEFAULT, MEM_DEFAULT, ROUTE_DEFAULT,
    )

    gcfg    = cfg.get("gaussian_opt", {}) or {}
    opt_com = work / f"{lig}_opt.com"

    work.mkdir(parents=True, exist_ok=True)

    try:
        write_com(
            str(ligand_pdb), str(opt_com),
            charge=charge, mult=mult,
            nproc=gcfg.get("nproc", NPROC_DEFAULT),
            mem=gcfg.get("mem",     MEM_DEFAULT),
            route=gcfg.get("route", ROUTE_DEFAULT),
        )
    except Exception as exc:
        return False, f"[{lig}] failed writing Gaussian opt input: {exc}"

    try:
        from senda.ligand_parameters.resp import write_resp_in, write_resp_qin
        write_resp_in(str(ligand_pdb), charge, resname, str(work / "resp.in"))
        write_resp_qin(str(ligand_pdb), str(work / "resp.qin"))
    except ImportError:
        return False, (
            f"[{lig}] RESP mode requires numpy (not available in this environment).\n"
            f"       Install with: conda install -c conda-forge numpy"
        )
    except Exception as exc:
        return False, f"[{lig}] failed writing RESP input files: {exc}"

    return True, (
        f"[{lig}] RESP Phase 1 done  →  {work}/\n"
        f"  Next steps:\n"
        f"    1. Submit {opt_com.name} to Gaussian on HPC\n"
        f"    2. Copy {lig}_opt.log back to {work}/\n"
        f"    3. Run: ligand_parameters-hf-input {work}/{lig}_opt.log\n"
        f"    4. Submit {lig}_hf.com to Gaussian on HPC\n"
        f"    5. Run: ligand_parameters-amber {work}/{lig}_hf.log"
    )


# ---------------------------------------------------------------------------
# Single-ligand dispatcher
# ---------------------------------------------------------------------------

def _run_single(lig: str, cwd: Path, cfg: dict, force: bool) -> tuple[bool, str]:
    from senda.ligand_parameters._utils import get_resname, resolve_charge

    ligands_dir = cwd / cfg.get("ligands_dir", _LIGANDS_DIR)
    output_dir  = cwd / cfg.get("output_dir",  _OUTPUT_DIR)
    logs_dir    = cwd / "logs"

    ligand_pdb = ligands_dir / f"{lig}.pdb"
    work       = output_dir / lig

    charge  = resolve_charge(ligand_pdb, cfg)
    resname = get_resname(ligand_pdb)
    method  = cfg.get("charge_method", _CHARGE_METHOD)
    mult    = cfg.get("multiplicity",  _MULTIPLICITY)
    ff      = (cfg.get("forcefield") or {}).get("ligand", _FORCEFIELD)

    # ── RESP: generate Gaussian inputs (Phase 1 only) ─────────────────────
    if method == "resp":
        opt_com = work / f"{lig}_opt.com"
        if not force and opt_com.exists():
            return True, f"[{lig}] RESP Phase 1 already done — skipping (use --force to redo)"
        return _run_resp_phase1(lig, cwd, cfg, work, ligand_pdb, charge, resname, mult)

    # ── BCC / gas: antechamber → parmchk2 → tleap ────────────────────────
    mol2       = work / f"{lig}.mol2"
    frcmod     = work / f"{lig}.frcmod"
    lib        = work / f"{lig}.lib"
    timing_log = work / "timing.log"
    local_log  = logs_dir / f"{lig}.log"

    if not force and mol2.exists() and frcmod.exists() and lib.exists():
        return True, f"[{lig}] already done — skipping (use --force to redo)"

    work.mkdir(parents=True, exist_ok=True)
    logs_dir.mkdir(parents=True, exist_ok=True)
    timing_log.unlink(missing_ok=True)

    tleap_in = _write_tleap_input(work, lig, resname, ff)

    t0 = t1 = time.monotonic()

    def _tick(label: str) -> None:
        nonlocal t1
        now = time.monotonic()
        ts  = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")
        with open(timing_log, "a") as fh:
            fh.write(f"{label}\t{int(now - t1)}\t{int(now - t0)}\t{ts}\n")
        t1 = now

    _tick("start")

    steps = [
        ("antechamber", [
            "antechamber",
            "-i", str(ligand_pdb), "-fi", "pdb",
            "-o", str(mol2),       "-fo", "mol2",
            "-c", method, "-s", "2",
            "-nc", str(charge), "-m", str(mult), "-rn", resname,
        ]),
        ("parmchk2", [
            "parmchk2", "-i", str(mol2), "-f", "mol2", "-o", str(frcmod),
        ]),
        ("tleap", [
            "tleap", "-f", str(tleap_in),
        ]),
    ]

    with open(local_log, "w") as fh:
        for step, cmd in steps:
            log.info("[%s] %s", lig, step)
            result = subprocess.run(cmd, cwd=str(work), stdout=fh, stderr=subprocess.STDOUT)
            _tick(step)
            if result.returncode != 0:
                tail = "\n".join(Path(local_log).read_text().splitlines()[-30:])
                return False, (
                    f"[{lig}] {step} failed (exit {result.returncode})\n"
                    f"Log: {local_log}\n--- last lines ---\n{tail}"
                )

    return True, f"[{lig}] done  →  {work}/"


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def run(cwd: Path, ligs: list[str], cfg: dict,
        jobs: int = 1, force: bool = False) -> None:
    ligands_dir = cwd / cfg.get("ligands_dir", _LIGANDS_DIR)
    output_dir  = cwd / cfg.get("output_dir",  _OUTPUT_DIR)
    method      = cfg.get("charge_method", _CHARGE_METHOD)

    print(f"Charge method : {method}")
    print(f"Ligands dir   : {ligands_dir}")
    print(f"Output dir    : {output_dir}")
    print(f"Ligands       : {ligs}")
    print(f"Jobs          : {jobs}")
    print()

    failed: list[str] = []

    if jobs == 1:
        for lig in ligs:
            ok, msg = _run_single(lig, cwd, cfg, force)
            print(msg)
            if not ok:
                failed.append(lig)
    else:
        with concurrent.futures.ThreadPoolExecutor(max_workers=jobs) as pool:
            futures = {pool.submit(_run_single, lig, cwd, cfg, force): lig
                       for lig in ligs}
            for fut in concurrent.futures.as_completed(futures):
                ok, msg = fut.result()
                print(msg)
                if not ok:
                    failed.append(futures[fut])

    if failed:
        log.error("Failed: %s", failed)
        sys.exit(1)

    if method == "resp":
        print(f"\nPhase 1 complete. Submit the Gaussian jobs above, then continue with "
              f"ligand_parameters-hf-input and ligand_parameters-amber.")
    else:
        print(f"\nAll done. Outputs in: {output_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main() -> None:
    import argparse

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    parser = argparse.ArgumentParser(
        description=(
            "Parameterize all ligands in ./ligands/ using AMBER/GAFF.\n"
            "Charge method is read from config.yaml (bcc or resp).\n"
            "Outputs go to ./ligands_libraries/<ligand>/."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", required=True, metavar="CONFIG",
        help="Path to the senda config.yaml")
    parser.add_argument("--ligands", nargs="+", default=None,
        help="Ligand IDs to process (default: all PDBs in ligands/)")
    parser.add_argument("--jobs", type=int, default=1, metavar="N",
        help="Parallel ligands (default: 1)")
    parser.add_argument("--force", action="store_true",
        help="Re-run even if outputs already exist")
    args = parser.parse_args()

    cwd = Path.cwd()
    cfg, inhibitors = _load_config(Path(args.config))

    ligands_dir = cwd / cfg.get("ligands_dir", _LIGANDS_DIR)
    if not ligands_dir.exists():
        sys.exit(
            f"Error: ligands directory not found: {ligands_dir}\n"
            f"Create a '{_LIGANDS_DIR}/' folder with .pdb files and run senda-param from there."
        )

    if args.ligands:
        ligs = args.ligands
    elif inhibitors:
        ligs = inhibitors
        log.info("Ligand list from config inhibitors: %s", ligs)
    else:
        ligs = _discover_ligands(ligands_dir)

    # Validate every ligand has a PDB before starting any work
    missing = [l for l in ligs if not (ligands_dir / f"{l}.pdb").exists()]
    if missing:
        sys.exit(
            f"Error: no PDB file found in {ligands_dir} for: {', '.join(missing)}"
        )

    run(cwd, ligs, cfg, jobs=args.jobs, force=args.force)

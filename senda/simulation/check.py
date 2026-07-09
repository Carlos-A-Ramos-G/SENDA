"""
senda.simulation.check
Report completion status and first failure point for all simulation replicas.

For each replica, stages are inspected in order (00_prep -> 01_min -> 02_heat
-> 03_equil -> 04_NVT).  The first stage that is incomplete or missing is
reported as the failure point.

Success criterion for each stage:
  00_prep  : structure.parm7 exists
  01_min   : last structure_min_N.out has the AMBER completion marker
  02_heat  : structure_heat.out has the AMBER completion marker
  03_equil : all 6 structure_equil_N.out files have the completion marker
  04_NVT   : structure_NVT_N.out files up to total_chunks are complete
"""
from __future__ import annotations

from pathlib import Path

_COMPLETION_MARKER = "Total wall time:"
_TAIL_BYTES        = 512   # how many bytes to read from the end of a .out file


def _completed(out: Path) -> bool:
    """Return True if an AMBER .out file ends with the completion marker."""
    if not out.exists() or out.stat().st_size == 0:
        return False
    with out.open("rb") as fh:
        fh.seek(max(0, out.stat().st_size - _TAIL_BYTES))
        return _COMPLETION_MARKER in fh.read().decode("utf-8", errors="replace")


def _tail(path: Path, n: int = 10) -> list[str]:
    """Return the last n non-empty lines of a text file."""
    if not path.exists():
        return ["(file not found)"]
    lines = [l for l in path.read_text(errors="replace").splitlines() if l.strip()]
    return lines[-n:] if lines else ["(empty file)"]


def _slurm_err(replica_dir: Path) -> Path | None:
    errs = sorted(replica_dir.glob("slurm-*.err"))
    return errs[-1] if errs else None


def _check_replica(replica_dir: Path, total_nvt_chunks: int) -> dict:
    """
    Returns a dict with keys:
      stage      : label of the current/failed stage
      status     : "complete" | "running" | "failed" | "not_started"
      progress   : human-readable progress string
      error_file : Path to the most informative file for diagnosis, or None
    """
    prep_dir  = replica_dir / "00_prep"
    min_dir   = replica_dir / "01_min"
    heat_dir  = replica_dir / "02_heat"
    equil_dir = replica_dir / "03_equil"
    nvt_dir   = replica_dir / "04_NVT"

    # --- 00_prep -------------------------------------------------------------
    if not (prep_dir / "structure.parm7").exists():
        return {
            "stage":      "00_prep",
            "status":     "failed",
            "progress":   "tleap",
            "error_file": _slurm_err(replica_dir),
        }

    # --- 01_min --------------------------------------------------------------
    min_outs = sorted(min_dir.glob("structure_min_[1-9]*.out"))
    if not min_outs:
        return {
            "stage":      "01_min",
            "status":     "not_started",
            "progress":   "0 cycles",
            "error_file": _slurm_err(replica_dir),
        }
    n_min_ok = sum(1 for f in min_outs if _completed(f))
    if not _completed(min_outs[-1]):
        return {
            "stage":      "01_min",
            "status":     "failed",
            "progress":   f"{n_min_ok}/{len(min_outs)} cycles",
            "error_file": min_outs[-1],
        }

    # --- 02_heat -------------------------------------------------------------
    heat_out = heat_dir / "structure_heat.out"
    if not heat_out.exists():
        return {
            "stage":      "02_heat",
            "status":     "not_started",
            "progress":   "-",
            "error_file": _slurm_err(replica_dir),
        }
    if not _completed(heat_out):
        return {
            "stage":      "02_heat",
            "status":     "failed",
            "progress":   "-",
            "error_file": heat_out,
        }

    # --- 03_equil ------------------------------------------------------------
    equil_outs = [equil_dir / f"structure_equil_{i}.out" for i in range(1, 7)]
    n_equil_ok = sum(1 for f in equil_outs if _completed(f))
    if n_equil_ok < 6:
        failing = next((f for f in equil_outs if not _completed(f)), None)
        status  = "failed" if (failing and failing.exists()) else "not_started"
        return {
            "stage":      "03_equil",
            "status":     status,
            "progress":   f"{n_equil_ok}/6",
            "error_file": failing if (failing and failing.exists()) else None,
        }

    # --- 04_NVT --------------------------------------------------------------
    nvt_outs = sorted(nvt_dir.glob("structure_NVT_[1-9]*.out"))
    n_nvt_ok = sum(1 for f in nvt_outs if _completed(f))
    total_label = f"/{total_nvt_chunks}" if total_nvt_chunks else ""

    if total_nvt_chunks and n_nvt_ok == total_nvt_chunks:
        return {
            "stage":      "04_NVT",
            "status":     "complete",
            "progress":   f"{n_nvt_ok}{total_label}",
            "error_file": None,
        }
    if nvt_outs and not _completed(nvt_outs[-1]):
        return {
            "stage":      "04_NVT",
            "status":     "failed",
            "progress":   f"{n_nvt_ok}{total_label}",
            "error_file": nvt_outs[-1],
        }
    return {
        "stage":      "04_NVT",
        "status":     "running",
        "progress":   f"{n_nvt_ok}{total_label}",
        "error_file": None,
    }


def check_all(
    inhibitors:      list[str],
    mutants:         list[str],
    n_replicas:      int,
    sim:             dict,
    simulations_dir: Path,
    verbose:         bool = False,
    only_failed:     bool = False,
) -> None:
    total_nvt_chunks = int((sim.get("production") or {}).get("total_chunks", 0))

    label_w  = max(
        len(f"{inh}/{mut}/replica_{n_replicas}")
        for inh in inhibitors for mut in mutants
    ) + 2
    stage_w  = 10
    status_w = 12
    progress_w = 12

    header = (f"  {'System':<{label_w}} {'Stage':<{stage_w}} "
              f"{'Status':<{status_w}} {'Progress':<{progress_w}}")
    print(header)
    print("  " + "-" * (len(header) - 2))

    counts = {"complete": 0, "running": 0, "failed": 0, "not_started": 0, "missing": 0}

    for inh in inhibitors:
        for mut in mutants:
            for rep in range(1, n_replicas + 1):
                label       = f"{inh}/{mut}/replica_{rep}"
                replica_dir = simulations_dir / inh / mut / f"replica_{rep}"

                if not replica_dir.exists():
                    counts["missing"] += 1
                    if not only_failed:
                        print(f"  {label:<{label_w}} {'--':<{stage_w}} "
                              f"{'missing':<{status_w}} --")
                    continue

                r       = _check_replica(replica_dir, total_nvt_chunks)
                status  = r["status"]
                stage   = r["stage"]
                progress = r["progress"]

                counts[status] = counts.get(status, 0) + 1

                is_bad = status in ("failed", "not_started")
                if only_failed and not is_bad:
                    continue

                flag = "* " if is_bad else "  "
                print(f"{flag}{label:<{label_w}} {stage:<{stage_w}} "
                      f"{status:<{status_w}} {progress:<{progress_w}}")

                if verbose and r["error_file"]:
                    err_path = r["error_file"]
                    print(f"    --> {err_path.relative_to(simulations_dir)}")
                    for line in _tail(err_path, 10):
                        print(f"      {line}")
                    print()

    print()
    parts = [f"{v} {k}" for k, v in counts.items() if v]
    print("  Summary: " + ", ".join(parts))

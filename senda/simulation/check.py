"""
senda.simulation.check
Report progress status for all simulation replicas.

For each replica, stages are inspected in order:
  00_prep -> 01_min -> 02_heat -> 03_equil -> 04_NVT

Output shows the last fully completed stage and the stage currently in progress
(if any).  A truncated .out file means the simulation may still be running --
it is never labelled as "failed".

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
_TAIL_BYTES        = 512


def _completed(out: Path) -> bool:
    if not out.exists() or out.stat().st_size == 0:
        return False
    with out.open("rb") as fh:
        fh.seek(max(0, out.stat().st_size - _TAIL_BYTES))
        return _COMPLETION_MARKER in fh.read().decode("utf-8", errors="replace")


def _parse_nstep(out: Path) -> int | None:
    """Return the last NSTEP value recorded in an AMBER .out file."""
    if not out.exists() or out.stat().st_size == 0:
        return None
    with out.open("rb") as fh:
        fh.seek(max(0, out.stat().st_size - 8192))
        text = fh.read().decode("utf-8", errors="replace")
    for line in reversed(text.splitlines()):
        if " NSTEP =" in line:
            try:
                return int(line.split("NSTEP =")[1].split()[0])
            except (IndexError, ValueError):
                pass
    return None


def _parse_nstlim(inp: Path) -> int | None:
    """Return nstlim from an AMBER .in input file."""
    if not inp.exists():
        return None
    for line in inp.read_text(errors="replace").splitlines():
        if "nstlim" in line.lower():
            try:
                return int(line.split("=")[1].strip().rstrip(",").split()[0])
            except (IndexError, ValueError):
                pass
    return None


def _tail(path: Path, n: int = 10) -> list[str]:
    if not path.exists():
        return ["(file not found)"]
    lines = [l for l in path.read_text(errors="replace").splitlines() if l.strip()]
    return lines[-n:] if lines else ["(empty file)"]


def _slurm_err(replica_dir: Path) -> Path | None:
    errs = sorted(replica_dir.glob("slurm-*.err"))
    return errs[-1] if errs else None


def _check_replica(replica_dir: Path, total_nvt_chunks: int) -> dict:
    """
    Returns:
      last_done      : label of last fully completed stage, or None
      last_done_prog : progress string for that stage
      in_progress    : label of the stage currently in progress, or None
      in_prog_prog   : progress string for in-progress stage
      active_file    : Path to the in-progress .out file (for -v), or None
    """
    prep_dir  = replica_dir / "00_prep"
    min_dir   = replica_dir / "01_min"
    heat_dir  = replica_dir / "02_heat"
    equil_dir = replica_dir / "03_equil"
    nvt_dir   = replica_dir / "04_NVT"

    last_done      = None
    last_done_prog = ""

    def _ret(in_progress=None, in_prog_prog="", active_file=None):
        return dict(
            last_done=last_done, last_done_prog=last_done_prog,
            in_progress=in_progress, in_prog_prog=in_prog_prog,
            active_file=active_file,
        )

    # 00_prep
    if not (prep_dir / "structure.parm7").exists():
        return _ret()
    last_done = "00_prep"

    # 01_min
    min_outs = sorted(min_dir.glob("structure_min_[1-9]*.out"))
    if not min_outs:
        return _ret()
    n_min_ok = sum(1 for f in min_outs if _completed(f))
    if _completed(min_outs[-1]):
        last_done = "01_min"
        last_done_prog = f"{len(min_outs)} stages"
    else:
        return _ret("01_min", f"{n_min_ok}/{len(min_outs)} stages", min_outs[-1])

    # 02_heat
    heat_out = heat_dir / "structure_heat.out"
    if not heat_out.exists():
        return _ret()
    if _completed(heat_out):
        last_done = "02_heat"
        nstlim = _parse_nstlim(heat_dir / "heat.in")
        last_done_prog = f"{nstlim} steps" if nstlim else "ok"
    else:
        nstep  = _parse_nstep(heat_out)
        nstlim = _parse_nstlim(heat_dir / "heat.in")
        if nstep is not None and nstlim:
            prog = f"{nstep}/{nstlim} steps"
        elif nstep is not None:
            prog = f"{nstep} steps"
        else:
            prog = ""
        return _ret("02_heat", prog, heat_out)

    # 03_equil
    equil_outs = [equil_dir / f"structure_equil_{i}.out" for i in range(1, 7)]
    n_equil_ok = sum(1 for f in equil_outs if _completed(f))
    if n_equil_ok == 6:
        last_done = "03_equil"
        last_done_prog = "6/6"
    elif any(f.exists() for f in equil_outs):
        active = next((f for f in equil_outs if f.exists() and not _completed(f)), None)
        return _ret("03_equil", f"{n_equil_ok}/6", active)
    else:
        return _ret()

    # 04_NVT
    nvt_outs = sorted(nvt_dir.glob("structure_NVT_[1-9]*.out"))
    n_nvt_ok = sum(1 for f in nvt_outs if _completed(f))
    total_label = f"/{total_nvt_chunks}" if total_nvt_chunks else ""

    is_done = (total_nvt_chunks and n_nvt_ok == total_nvt_chunks) or (
        not total_nvt_chunks and nvt_outs and _completed(nvt_outs[-1])
    )
    if is_done:
        last_done = "04_NVT"
        last_done_prog = f"{n_nvt_ok}{total_label}"
        return _ret()
    if nvt_outs:
        return _ret("04_NVT", f"{n_nvt_ok}{total_label}", nvt_outs[-1])

    return _ret()


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

    label_w = max(
        len(f"{inh}/{mut}/replica_{n_replicas}")
        for inh in inhibitors for mut in mutants
    ) + 2
    comp_w = 22   # "Completed" column width (stage + progress)

    header = f"  {'System':<{label_w}} {'Completed':<{comp_w}} In progress"
    print(header)
    print("  " + "-" * (len(header) - 2))

    counts = {"complete": 0, "in_progress": 0, "not_started": 0, "missing": 0}

    for inh in inhibitors:
        for mut in mutants:
            for rep in range(1, n_replicas + 1):
                label       = f"{inh}/{mut}/replica_{rep}"
                replica_dir = simulations_dir / inh / mut / f"replica_{rep}"

                if not replica_dir.exists():
                    counts["missing"] += 1
                    if not only_failed:
                        print(f"  {label:<{label_w}} {'--':<{comp_w}} --")
                    continue

                r = _check_replica(replica_dir, total_nvt_chunks)
                last_done   = r["last_done"]
                in_progress = r["in_progress"]

                is_complete = (last_done == "04_NVT" and in_progress is None)
                if is_complete:
                    counts["complete"] += 1
                elif in_progress:
                    counts["in_progress"] += 1
                elif last_done:
                    counts["not_started"] += 1   # paused between stages
                else:
                    counts["not_started"] += 1

                if only_failed and is_complete:
                    continue

                comp_str = (
                    f"{last_done}  {r['last_done_prog']}".rstrip()
                    if last_done else "--"
                )
                prog_str = (
                    f"{in_progress}  {r['in_prog_prog']}".rstrip()
                    if in_progress else "--"
                )
                print(f"  {label:<{label_w}} {comp_str:<{comp_w}} {prog_str}")

                if verbose and r["active_file"]:
                    af = r["active_file"]
                    print(f"    --> {af.relative_to(simulations_dir)}")
                    for line in _tail(af, 10):
                        print(f"      {line}")
                    print()

    print()
    parts = [f"{v} {k}" for k, v in counts.items() if v]
    print("  Summary: " + ", ".join(parts))

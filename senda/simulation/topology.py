"""
senda.simulation.topology
Print a summary table of atoms, box dimensions, and ion counts for each
system, read from replica_1/00_prep/structure.parm7.
"""
from __future__ import annotations

from pathlib import Path

_WATER_NAMES = {"WAT", "HOH", "TIP3", "TP3", "T3P"}
_NA_NAMES    = {"Na+", "NA",  "SOD",  "NAP"}
_CL_NAMES    = {"Cl-", "CL",  "CLA",  "CLM"}


def _read_parm7(parm7: Path) -> dict | None:
    try:
        import parmed as pmd
    except ImportError:
        raise SystemExit("ParmEd is required: pip install parmed  (or load AmberTools)")

    try:
        top = pmd.load_file(str(parm7))
    except Exception:
        return None

    box      = top.box   # [a, b, c, alpha, beta, gamma] in Angstroms/degrees, or None
    ion_solv = _WATER_NAMES | _NA_NAMES | _CL_NAMES
    names    = [r.name for r in top.residues]

    solute_q = round(sum(
        a.charge for r in top.residues if r.name not in ion_solv for a in r.atoms
    ))

    return {
        "atoms":    len(top.atoms),
        "box_x":    float(box[0]) if box is not None else None,
        "box_y":    float(box[1]) if box is not None else None,
        "box_z":    float(box[2]) if box is not None else None,
        "waters":   sum(1 for n in names if n in _WATER_NAMES),
        "na":       sum(1 for n in names if n in _NA_NAMES),
        "cl":       sum(1 for n in names if n in _CL_NAMES),
        "solute_q": solute_q,
    }


def topology_info(
    inhibitors:      list[str],
    mutants:         list[str],
    simulations_dir: Path,
) -> None:
    rows = []
    for inh in inhibitors:
        for mut in mutants:
            parm7 = simulations_dir / inh / mut / "replica_1" / "00_prep" / "structure.parm7"
            label = f"{inh}/{mut}"
            rows.append((label, _read_parm7(parm7) if parm7.exists() else None))

    sys_w = max(max(len(r[0]) for r in rows), len("System"))

    header = (
        f"  {'System':<{sys_w}}  {'Atoms':>8}  "
        f"{'Box X (A)':>9}  {'Box Y (A)':>9}  {'Box Z (A)':>9}  "
        f"{'Waters':>7}  {'Na+':>4}  {'Cl-':>4}  {'Solute Q':>8}"
    )
    print(header)
    print("  " + "-" * (len(header) - 2))

    for label, info in rows:
        if info is None:
            dash = "--"
            print(
                f"  {label:<{sys_w}}  {dash:>8}  "
                f"{dash:>9}  {dash:>9}  {dash:>9}  "
                f"{dash:>7}  {dash:>4}  {dash:>4}  {dash:>8}"
            )
            continue

        box_x = f"{info['box_x']:.2f}" if info["box_x"] is not None else "--"
        box_y = f"{info['box_y']:.2f}" if info["box_y"] is not None else "--"
        box_z = f"{info['box_z']:.2f}" if info["box_z"] is not None else "--"

        print(
            f"  {label:<{sys_w}}  {info['atoms']:>8,}  "
            f"{box_x:>9}  {box_y:>9}  {box_z:>9}  "
            f"{info['waters']:>7,}  {info['na']:>4}  {info['cl']:>4}  {info['solute_q']:>+8}"
        )

"""
senda.michaelis_complex.disulfides
Detect candidate Cys SG-SG disulfide bonds by distance, for the
senda-complex --report-disulfides preliminary run.
"""

from __future__ import annotations

import numpy as np

from .clean import _chain, _coords, _resname, _resnum

_BONDED_CUTOFF    = 2.5   # Angstrom -- typical S-S bond length is ~2.05 A
_BORDERLINE_CUTOFF = 3.0  # Angstrom -- flagged for manual review, not auto-suggested


def _atom_name(line: str) -> str:
    return line[12:16].strip()


def find_disulfide_candidates(
    lines:          list[str],
    chains:         frozenset[str],
    bonded_cutoff:  float = _BONDED_CUTOFF,
    check_cutoff:   float = _BORDERLINE_CUTOFF,
) -> list[dict]:
    """
    Return candidate Cys SG-SG disulfide pairs found in *lines*.

    Each result is a dict with keys chain1, resnum1, chain2, resnum2, distance,
    band ("bonded" if distance <= bonded_cutoff, else "borderline"). Pairs
    farther apart than check_cutoff are not reported at all.

    Only the primary altloc conformer (blank or "A") is considered for each
    (chain, resnum), so a split-occupancy Cys contributes a single SG position.
    """
    sg_atoms: dict[tuple[str, int], tuple[float, float, float]] = {}

    for line in lines:
        if not line.startswith("ATOM  "):
            continue
        if _resname(line) not in ("CYS", "CYX"):
            continue
        if _atom_name(line) != "SG":
            continue
        altloc = line[16]
        if altloc not in (" ", "A"):
            continue
        ch = _chain(line)
        if ch not in chains:
            continue
        rn = _resnum(line)
        if rn is None:
            continue
        xyz = _coords(line)
        if xyz is None:
            continue
        sg_atoms.setdefault((ch, rn), xyz)

    keys = list(sg_atoms.keys())
    if len(keys) < 2:
        return []

    coords = np.asarray([sg_atoms[k] for k in keys])
    dists  = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

    candidates: list[dict] = []
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            d = float(dists[i, j])
            if d > check_cutoff:
                continue
            ch1, rn1 = keys[i]
            ch2, rn2 = keys[j]
            candidates.append({
                "chain1":   ch1,
                "resnum1":  rn1,
                "chain2":   ch2,
                "resnum2":  rn2,
                "distance": d,
                "band":     "bonded" if d <= bonded_cutoff else "borderline",
            })

    candidates.sort(key=lambda c: c["distance"])
    return candidates

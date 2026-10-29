"""
senda.michaelis_complex.align
Cα superposition of a mobile PDB onto a reference PDB.

Uses Biopython SVDSuperimposer for the Kabsch rotation (rotation matrix R and
translation vector t).  The transform is then applied directly to PDB coordinate
columns (cols 30-54) of every ATOM/HETATM record.  ANISOU records are left
unchanged (they are stripped by the simulation launcher anyway).
"""

from __future__ import annotations

from senda.config import AMINO_ACIDS


# ---------------------------------------------------------------------------
# Cα extraction
# ---------------------------------------------------------------------------

def _parse_ca(
    lines:         list[str],
    chains:        tuple[str, ...]        = ("A",),
    resnum_range:  tuple[int, int] | None = None,
) -> tuple[list[tuple[str, int]], list[list[float]]]:
    """
    Return (keys, coords) for the CA atoms of amino-acid residues in *chains*
    (optionally only residue numbers within *resnum_range*, inclusive).
    key = (chain_id, resnum); a blank chain ID counts as A.

    Only amino-acid residues count, so a calcium ion or a ligand atom named
    CA never enters the fit. Aligning on a single chain (the default chain)
    keeps multimeric crystal structures from being misaligned by another
    chain whose relative orientation differs from the reference.
    """
    keys:   list[tuple[str, int]] = []
    coords: list[list[float]]     = []

    for line in lines:
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        if line[17:20].strip() not in AMINO_ACIDS:
            continue
        chain = line[21].strip() or "A"
        if chain not in chains:
            continue
        try:
            resnum = int(line[22:26])
            x      = float(line[30:38])
            y      = float(line[38:46])
            z      = float(line[46:54])
        except ValueError:
            continue
        if resnum_range is not None and not (resnum_range[0] <= resnum <= resnum_range[1]):
            continue
        keys.append((chain, resnum))
        coords.append([x, y, z])

    return keys, coords


def _matched_arrays(
    mob_keys: list[tuple[str, int]], mob_coords: list[list[float]],
    ref_keys: list[tuple[str, int]], ref_coords: list[list[float]],
) -> tuple[list[list[float]], list[list[float]]]:
    """Return paired (mobile, reference) Cα coordinate lists for shared keys."""
    ref_map = {k: c for k, c in zip(ref_keys, ref_coords)}
    mob_xyz: list[list[float]] = []
    ref_xyz: list[list[float]] = []

    for k, c in zip(mob_keys, mob_coords):
        if k in ref_map:
            mob_xyz.append(c)
            ref_xyz.append(ref_map[k])

    if not mob_xyz:
        raise ValueError(
            "No shared Cα atoms found between mobile and reference structures. "
            "Check michaelis_complex.align_chains / align_residues: both files "
            "need amino-acid CA atoms in those chains and residue numbers."
        )

    return mob_xyz, ref_xyz


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def compute_superposition(
    mobile_lines: list[str],
    ref_lines:    list[str],
    chains:       tuple[str, ...]        = ("A",),
    resnum_range: tuple[int, int] | None = None,
) -> tuple:
    """
    Compute (R, t, rms) such that  R @ mobile_coord + t ≈ ref_coord.

    Uses Biopython's SVDSuperimposer (Kabsch algorithm).
    Requires: biopython + numpy.

    Returns
    -------
    rot  : (3, 3) rotation matrix  (numpy array)
    tran : (3,)   translation vector (numpy array)
    rms  : float  Cα RMSD after superposition
    """
    try:
        import numpy as np
        from Bio.SVDSuperimposer import SVDSuperimposer
    except ImportError as exc:
        raise ImportError(
            "biopython and numpy are required for Michaelis complex preparation.\n"
            "Install with:  conda install -c conda-forge biopython numpy"
        ) from exc

    mob_keys, mob_coords = _parse_ca(mobile_lines, chains, resnum_range)
    ref_keys, ref_coords = _parse_ca(ref_lines, chains, resnum_range)

    mob_xyz, ref_xyz = _matched_arrays(mob_keys, mob_coords, ref_keys, ref_coords)

    mob_arr = np.array(mob_xyz, dtype=float)
    ref_arr = np.array(ref_xyz, dtype=float)

    sup = SVDSuperimposer()
    sup.set(ref_arr, mob_arr)   # set(reference, mobile)
    sup.run()

    rot, tran = sup.get_rotran()
    rms       = sup.get_rms()

    return rot, tran, rms


def apply_transform(lines: list[str], rot, tran) -> list[str]:
    """
    Apply rotation + translation to x/y/z columns of every ATOM/HETATM record.
    ANISOU records are passed through unchanged.
    rot and tran are numpy arrays returned by compute_superposition().
    """
    import numpy as np

    out: list[str] = []

    for line in lines:
        if line.startswith(("ATOM  ", "HETATM")):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                out.append(line)
                continue
            nx, ny, nz = np.array([x, y, z]) @ rot + tran
            line = (
                line[:30]
                + f"{nx:8.3f}{ny:8.3f}{nz:8.3f}"
                + line[54:]
            )
        out.append(line)

    return out

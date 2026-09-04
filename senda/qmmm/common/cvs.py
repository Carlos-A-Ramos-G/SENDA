"""
senda.qmmm.common.cvs

CVs file writer and guess file interpolation / IO.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np


# AMBER COLVAR_type names per CV type
_AMBER_CV_TYPE = {
    "distance": "BOND",
    "angle":    "ANGLE",
    "dihedral": "TORSION",
}

# AMBER r-bounds for each CV type (before/after the target window)
TAIL_DIST  = 999.0  # Angstrom tails for distance/angle
TAIL_DIHED = 180.0  # degrees added to dihedral target for tail


# ---------------------------------------------------------------------------
# CV restraint (&rst) blocks
# ---------------------------------------------------------------------------

def build_rst_block(
    indices:        list[int],
    target:         float,
    cv_type:        str,
    force_constant: float,
) -> str:
    """Return a single AMBER &rst block restraining one CV to target."""
    iat = ", ".join(str(i) for i in indices) + ","
    fc  = force_constant

    if cv_type == "distance":
        r1 = 0.0
        r2 = target
        r3 = target
        r4 = TAIL_DIST
    elif cv_type == "angle":
        r1 = max(0.0, target - 180.0)
        r2 = target
        r3 = target
        r4 = min(360.0, target + 180.0)
    elif cv_type == "dihedral":
        r1 = target - TAIL_DIHED
        r2 = target
        r3 = target
        r4 = target + TAIL_DIHED
    else:
        r1 = 0.0
        r2 = target
        r3 = target
        r4 = TAIL_DIST

    return (
        f"&rst\n"
        f" iat={iat}\n"
        f" r1={r1:.4f}, r2={r2:.4f}, r3={r3:.4f}, r4={r4:.4f},\n"
        f" rk2={fc:.2f}, rk3={fc:.2f},\n"
        f"/\n"
    )


# ---------------------------------------------------------------------------
# CVs file
# ---------------------------------------------------------------------------

def write_cvs_file(
    cv_specs:          list,
    cv_indices_per_cv: list[list[int]],
    out_path:          Path,
) -> None:
    """Write the AMBER-format CVs file used by sander's string method."""
    n_cvs = len(cv_specs)
    lines = [f"{n_cvs}\n\n"]
    for cv, indices in zip(cv_specs, cv_indices_per_cv):
        cv_type    = cv.get("type", "distance").lower()
        amber_type = _AMBER_CV_TYPE.get(cv_type, "BOND")
        atom_str   = ", ".join(str(i) for i in indices)
        lines.append("$COLVAR\n")
        lines.append(f'COLVAR_type = "{amber_type}"\n')
        lines.append(f"atoms = {atom_str}\n")
        lines.append("$END\n\n")
    out_path.write_text("".join(lines))


# ---------------------------------------------------------------------------
# Guess file IO
# ---------------------------------------------------------------------------

def load_guess(guess_path: Path) -> np.ndarray:
    """
    Load a guess file as a plain 2-D numpy array (nodes x n_cvs).
    Accepts files with or without a header line; strips comment lines.
    """
    rows = []
    with open(guess_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                values = [float(v) for v in line.split()]
                rows.append(values)
            except ValueError:
                continue  # skip non-numeric lines (e.g. header)

    if not rows:
        raise ValueError(f"Guess file {guess_path} contains no numeric data")

    # Remove the header row if present (AMBER string header: N_nodes  N_cvs  energy)
    # Heuristic: if the first row has fewer columns than the second, it is a header
    if len(rows) > 1 and len(rows[0]) < len(rows[1]):
        rows = rows[1:]

    arr = np.array(rows, dtype=float)
    if arr.ndim != 2:
        raise ValueError(f"Guess file {guess_path}: expected 2-D data, got shape {arr.shape}")
    return arr


def interpolate_guess(data: np.ndarray, n_nodes: int) -> np.ndarray:
    """
    Linearly interpolate a (M x n_cvs) guess to (n_nodes x n_cvs) using
    arc-length parametrisation.  Returns data unchanged if M == n_nodes.
    """
    M, n_cvs = data.shape
    if M == n_nodes:
        return data
    if M < 2:
        raise ValueError("Guess must have at least 2 rows to interpolate")

    # Arc-length along the path
    diffs      = np.diff(data, axis=0)                          # (M-1, n_cvs)
    seg_len    = np.sqrt((diffs ** 2).sum(axis=1))              # (M-1,)
    arc        = np.concatenate([[0.0], np.cumsum(seg_len)])    # (M,)
    total      = arc[-1]

    new_t      = np.linspace(0.0, total, n_nodes)
    new_data   = np.zeros((n_nodes, n_cvs))
    for cv_i in range(n_cvs):
        new_data[:, cv_i] = np.interp(new_t, arc, data[:, cv_i])

    return new_data


def write_scan_guess(data: np.ndarray, out_path: Path) -> None:
    """Write plain N_nodes x N_cvs matrix (no header) for the scan step."""
    lines = []
    for row in data:
        lines.append("  ".join(f"{v:.6f}" for v in row))
    out_path.write_text("\n".join(lines) + "\n")


def write_string_guess(data: np.ndarray, out_path: Path) -> None:
    """Write guess with AMBER string header: bare N_nodes on its own line."""
    n_nodes, n_cvs = data.shape
    lines = [str(n_nodes)]
    for row in data:
        lines.append("  ".join(f"{v:.6f}" for v in row))
    out_path.write_text("\n".join(lines) + "\n")

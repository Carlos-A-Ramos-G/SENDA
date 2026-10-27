"""
senda.simulation.openmm_engine.restraints
AMBER-equivalent restraints as OpenMM forces.

Positional restraints (AMBER ntr=1)
-----------------------------------
AMBER:  E = restraint_wt * |r - r_ref|^2   (kcal/mol/A^2, *no* factor 1/2)
OpenMM: CustomExternalForce "k_pos*periodicdistance(x,y,z,x0,y0,z0)^2"
        k_pos [kJ/mol/nm^2] = restraint_wt * KCAL_A2_TO_KJ_NM2

NMR restraints (AMBER nmropt=1 + DISANG &rst)
---------------------------------------------
Flat-bottomed well with linear tails, identical to AMBER's form:

  r <  r1        : rk2*(r1-r2)^2 + 2*rk2*(r1-r2)*(r-r1)   (linear)
  r1 <= r < r2   : rk2*(r-r2)^2
  r2 <= r <= r3  : 0
  r3 <  r <= r4  : rk3*(r-r3)^2
  r  >  r4       : rk3*(r4-r3)^2 + 2*rk3*(r4-r3)*(r-r4)   (linear)

2 atoms -> distance (A, kcal/mol/A^2); 3 atoms -> angle, 4 atoms -> dihedral
(r in degrees, rk in kcal/mol/rad^2). Dihedrals are translated by 2*pi so the
value lies within pi of the midpoint of r2..r3, as AMBER does.

Indices: AMBER (and the DISANG file) are 1-based; everything returned here is
0-based (OpenMM/ASH convention). Conversion happens only in parse_disang.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

import openmm as mm

KCAL_TO_KJ       = 4.184
KCAL_A2_TO_KJ_NM2 = KCAL_TO_KJ * 100.0   # kcal/mol/A^2 -> kJ/mol/nm^2

_SUPPORTED_RST_KEYS = {"iat", "r1", "r2", "r3", "r4", "rk2", "rk3"}

_FLAT_BOTTOM = (
    "select(step(r1-r), rk2*(r1-r2)^2 + 2*rk2*(r1-r2)*(r-r1),"
    " select(step(r2-r), rk2*(r-r2)^2,"
    " select(step(r3-r), 0,"
    " select(step(r4-r), rk3*(r-r3)^2,"
    " rk3*(r4-r3)^2 + 2*rk3*(r4-r3)*(r-r4)))))"
)

_PARAMS = ("r1", "r2", "r3", "r4", "rk2", "rk3")


# ---------------------------------------------------------------------------
# Positional restraints
# ---------------------------------------------------------------------------

def positional_force(indices, ref_positions_nm, weight_kcal: float) -> mm.CustomExternalForce:
    """Harmonic positional restraint on *indices* (0-based) toward *ref_positions_nm*.

    *ref_positions_nm* is indexable by atom index and yields (x, y, z) in nm
    (e.g. an (N, 3) array). The force constant is the global parameter
    ``k_pos`` so a restraint schedule only needs ``context.setParameter``.
    """
    force = mm.CustomExternalForce("k_pos*periodicdistance(x,y,z,x0,y0,z0)^2")
    force.addGlobalParameter("k_pos", weight_kcal * KCAL_A2_TO_KJ_NM2)
    for p in ("x0", "y0", "z0"):
        force.addPerParticleParameter(p)
    for i in indices:
        x, y, z = (float(c) for c in ref_positions_nm[i])
        force.addParticle(int(i), [x, y, z])
    force.setName("PositionalRestraint")
    return force


def read_mask_atoms(path: Path) -> list[int]:
    """0-based atom indices from a cpptraj ``mask ... maskout`` file.

    Columns: #Frame AtomNum Atom ResNum Res MolNum (AtomNum is 1-based).
    """
    indices = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        indices.append(int(line.split()[1]) - 1)
    if not indices:
        raise ValueError(f"no atoms selected in {path} -- check the restraint mask")
    return indices


# ---------------------------------------------------------------------------
# NMR restraints
# ---------------------------------------------------------------------------

def parse_disang(path: Path) -> list[dict]:
    """Parse an AMBER DISANG file made of ``&rst ... /`` blocks.

    Returns dicts with 0-based ``atoms`` and r1..r4, rk2, rk3 in AMBER units.
    Only the iat/r1-r4/rk2/rk3 form is supported (what restrainer.py writes);
    anything else (igr*, r0/k0, nstep*, ...) raises instead of being ignored.
    """
    text = "\n".join(
        line for line in Path(path).read_text().splitlines()
        if not line.lstrip().startswith(("#", "!"))
    )
    blocks = re.findall(r"&rst\b(.*?)(?:^\s*/|/\s*$|&end)", text,
                        flags=re.DOTALL | re.IGNORECASE | re.MULTILINE)
    restraints = []
    for body in blocks:
        keys = list(re.finditer(r"([A-Za-z]\w*)\s*=", body))
        fields = {}
        for k, nxt in zip(keys, keys[1:] + [None]):
            raw = body[k.end(): nxt.start() if nxt else len(body)]
            vals = [v for v in re.split(r"[,\s]+", raw.strip()) if v]
            fields[k.group(1).lower()] = vals
        unknown = set(fields) - _SUPPORTED_RST_KEYS
        if unknown:
            raise ValueError(
                f"{path}: unsupported &rst keys {sorted(unknown)} for the OpenMM "
                f"engine (supported: {sorted(_SUPPORTED_RST_KEYS)})"
            )
        missing = _SUPPORTED_RST_KEYS - set(fields)
        if missing:
            raise ValueError(f"{path}: &rst block missing {sorted(missing)}")
        atoms = [int(v) - 1 for v in fields["iat"] if int(v) > 0]
        if len(atoms) not in (2, 3, 4):
            raise ValueError(f"{path}: &rst needs 2-4 atoms, got iat={fields['iat']}")
        r = {"atoms": atoms}
        for p in _PARAMS:
            r[p] = float(fields[p][0])
        restraints.append(r)
    return restraints


def _to_openmm_units(r: dict) -> list[float]:
    """AMBER units -> OpenMM units for one restraint, ordered as _PARAMS."""
    if len(r["atoms"]) == 2:
        return [r[p] * 0.1 for p in ("r1", "r2", "r3", "r4")] + \
               [r[p] * KCAL_A2_TO_KJ_NM2 for p in ("rk2", "rk3")]
    return [math.radians(r[p]) for p in ("r1", "r2", "r3", "r4")] + \
           [r[p] * KCAL_TO_KJ for p in ("rk2", "rk3")]


def nmr_forces(restraints: list[dict]) -> list[mm.Force]:
    """One CustomCompoundBondForce per restraint kind present in *restraints*."""
    kinds = {
        2: ("NMRDistance", "r = distance(p1, p2)"),
        3: ("NMRAngle",    "r = angle(p1, p2, p3)"),
        # translate theta by 2*pi so it lies within pi of the r2..r3 midpoint
        4: ("NMRDihedral", "r = theta - 6.283185307179586*floor((theta - mid + 3.141592653589793)"
                           "/6.283185307179586); mid = 0.5*(r2 + r3);"
                           " theta = dihedral(p1, p2, p3, p4)"),
    }
    forces = []
    for n_atoms, (name, definition) in kinds.items():
        group = [r for r in restraints if len(r["atoms"]) == n_atoms]
        if not group:
            continue
        force = mm.CustomCompoundBondForce(n_atoms, f"{_FLAT_BOTTOM}; {definition}")
        for p in _PARAMS:
            force.addPerBondParameter(p)
        for r in group:
            force.addBond(r["atoms"], _to_openmm_units(r))
        force.setUsesPeriodicBoundaryConditions(True)
        force.setName(name)
        forces.append(force)
    return forces

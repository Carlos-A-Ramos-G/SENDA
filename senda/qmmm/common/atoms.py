"""
senda.qmmm.common.atoms

CV atom resolution, QM region auto-selection, and H10 hydrogen identification.

Atom specs use the same syntax as analysis.reactive_distances:
  {sequence: N, name: atomname}         -- protein residue by PDB resnum
  {substrate_sequence: N, name: aname}  -- ligand/substrate by AMBER resid proximity
  {nearest_water_to: <ref_spec>, name: atomname}  -- WAT closest to ref in rep. frame
"""
from __future__ import annotations

from collections import defaultdict
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Standard protein residue names (used to detect backbone atoms)
# ---------------------------------------------------------------------------

_PROTEIN_RESNAMES = frozenset({
    "ALA", "ARG", "ASN", "ASP", "CYS", "CYM", "CYX",
    "GLN", "GLU", "GLH", "GLY", "HIS", "HID", "HIE", "HIP",
    "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
    "ASH", "ACE", "NME",
})

_SOLVENT_RESNAMES = frozenset({"WAT", "HOH", "Na+", "Cl-", "K+", "IP", "SOD", "CLA"})


# ---------------------------------------------------------------------------
# CV atom resolution
# ---------------------------------------------------------------------------

def resolve_cv_atoms_per_cv(
    cv_specs:            list,
    chain_map:           dict,
    top_atoms:           list,
    top_residues:        list,
    substrate_chain_map: dict,
    chain:               str,
    rst7_coords:         np.ndarray,
) -> list[list[int]]:
    """
    Resolve all CV atom specs to 1-based AMBER atom indices.

    Returns a list-of-lists: one inner list per CV, containing the resolved
    atom indices for that CV in the given chain.
    """
    result = []
    for cv in cv_specs:
        cv_indices = [
            _resolve_single(
                spec, chain_map, top_atoms, top_residues,
                substrate_chain_map, chain, rst7_coords,
            )
            for spec in cv["atoms"]
        ]
        result.append(cv_indices)
    return result


def resolve_cv_atoms_flat(
    cv_specs:            list,
    chain_map:           dict,
    top_atoms:           list,
    top_residues:        list,
    substrate_chain_map: dict,
    chain:               str,
    rst7_coords:         np.ndarray,
) -> list[int]:
    """Flat list of all resolved 1-based atom indices across all CVs."""
    indices = []
    for per_cv in resolve_cv_atoms_per_cv(
        cv_specs, chain_map, top_atoms, top_residues,
        substrate_chain_map, chain, rst7_coords,
    ):
        indices.extend(per_cv)
    return indices


def _resolve_single(
    spec:                dict,
    chain_map:           dict,
    top_atoms:           list,
    top_residues:        list,
    substrate_chain_map: dict,
    chain:               str,
    rst7_coords:         np.ndarray,
) -> int:
    if "nearest_water_to" in spec:
        return _resolve_nearest_water(
            spec["nearest_water_to"], spec["name"],
            chain_map, top_atoms, top_residues,
            substrate_chain_map, chain, rst7_coords,
        )
    from senda.analysis.distances import _resolve_atom_index
    return _resolve_atom_index(chain_map, top_atoms, spec, chain, substrate_chain_map)


def _resolve_nearest_water(
    ref_spec:            dict,
    water_atom_name:     str,
    chain_map:           dict,
    top_atoms:           list,
    top_residues:        list,
    substrate_chain_map: dict,
    chain:               str,
    rst7_coords:         np.ndarray,
) -> int:
    """Return 1-based index of water_atom_name in the WAT residue whose O is
    closest to ref_spec in rst7_coords."""
    from senda.analysis.distances import _find_nearest_water_resid
    best_resid = _find_nearest_water_resid(
        ref_spec, chain_map, top_atoms, top_residues,
        substrate_chain_map, chain, rst7_coords,
    )

    for atom in top_atoms:
        if atom.resid == best_resid and atom.name == water_atom_name:
            return atom.index + 1

    available = [a.name for a in top_atoms if a.resid == best_resid]
    raise ValueError(
        f"Atom {water_atom_name!r} not found in nearest WAT residue "
        f"(resid {best_resid + 1}). Available: {available}"
    )


def resolve_extra_restraints(
    restraint_specs:     list,
    chain_map:           dict,
    top_atoms:           list,
    top_residues:        list,
    substrate_chain_map: dict,
    chain:               str,
    rst7_coords:         np.ndarray,
) -> list[dict]:
    """
    Resolve each restraint's 'atoms' entries to 1-based AMBER indices.

    Raw integers pass through unchanged (already an atom index).
    {sequence: ...} / {substrate_sequence: ...} / {nearest_water_to: ...}
    dict specs are resolved via the same machinery CV atoms use.
    """
    resolved = []
    for r in restraint_specs:
        new_r = dict(r)
        new_r["atoms"] = [
            a if isinstance(a, int) else _resolve_single(
                a, chain_map, top_atoms, top_residues,
                substrate_chain_map, chain, rst7_coords,
            )
            for a in r["atoms"]
        ]
        resolved.append(new_r)
    return resolved


# ---------------------------------------------------------------------------
# QM region auto-selection
# ---------------------------------------------------------------------------

_QMMASK_WATER_PLACEHOLDER = "__NEAREST_WATER__"


def _canon_water_ref(ref_spec) -> frozenset:
    """Normalize a nearest_water_to ref (single spec or list of specs) to an
    order-independent, hashable form so two refs can be compared for
    equality regardless of list ordering."""
    specs = ref_spec if isinstance(ref_spec, list) else [ref_spec]
    return frozenset(tuple(sorted(spec.items())) for spec in specs)


def _cv_water_refs(cv_specs: list) -> list:
    """Collect every nearest_water_to ref used by any CV atom."""
    refs = []
    for cv in cv_specs:
        for atom_spec in cv.get("atoms", []):
            if isinstance(atom_spec, dict) and "nearest_water_to" in atom_spec:
                refs.append(atom_spec["nearest_water_to"])
    return refs


def resolve_qm_region(
    inh_cfg:             dict,
    cv_indices_1based:   list[int],
    parm_path:           Path,
    chain_map:           dict | None = None,
    top_atoms:            list | None = None,
    top_residues:         list | None = None,
    substrate_chain_map: dict | None = None,
    chain:               str  | None = None,
    rst7_coords:          np.ndarray | None = None,
    cv_specs:             list | None = None,
) -> tuple[str, int]:
    """
    Return (qmmask, qmcharge).

    If 'qmmask' is present in inh_cfg: use it and require 'qmcharge'. If the
    mask string contains the '__NEAREST_WATER__' placeholder (write it as
    ':__NEAREST_WATER__' -- the leading ':' is part of the mask, not the
    placeholder), it is replaced with the residue number of the WAT residue
    nearest 'qmwater_neighbor' -- resolved fresh from rst7_coords, same as a
    CV's nearest_water_to, so it stays correct across mutants and
    re-selected representative frames. Requires the chain_map/top_atoms/
    top_residues/substrate_chain_map/chain/rst7_coords arguments in that
    case.

    'qmwater_neighbor' is optional: if every CV's nearest_water_to agrees on
    one water, that's used automatically -- this guarantees the QM region's
    water is the SAME one the CVs (and H10's mass patch) actually use,
    rather than requiring it to be hand-declared and kept in sync
    separately, which can silently drift out of sync (the QM region ending
    up with a different water than the one driving the reaction coordinate).
    If the CVs reference more than one distinct water, or none at all,
    'qmwater_neighbor' must be set explicitly; if it IS set, it must match
    one of the CVs' nearest_water_to refs (when any exist) or resolution is
    rejected with an error, rather than silently allowing the mismatch.

    Otherwise: auto-select from CV atom seed residues via bond-graph expansion.
    """
    if "qmmask" in inh_cfg:
        if "qmcharge" not in inh_cfg:
            raise ValueError(
                "'qmmask' is set manually but 'qmcharge' is not provided. "
                "When overriding qmmask, qmcharge must be declared explicitly."
            )
        qmmask = inh_cfg["qmmask"]
        if _QMMASK_WATER_PLACEHOLDER in qmmask:
            from senda.analysis.distances import _find_nearest_water_resid

            cv_refs = _cv_water_refs(cv_specs or [])
            distinct_cv_refs = {_canon_water_ref(r): r for r in cv_refs}

            water_ref = inh_cfg.get("qmwater_neighbor")
            if water_ref is None:
                if not distinct_cv_refs:
                    raise ValueError(
                        f"qmmask contains {_QMMASK_WATER_PLACEHOLDER!r} but "
                        "'qmwater_neighbor' is not set and no CV uses "
                        "nearest_water_to to infer it from."
                    )
                if len(distinct_cv_refs) > 1:
                    raise ValueError(
                        f"qmmask contains {_QMMASK_WATER_PLACEHOLDER!r} but "
                        "'qmwater_neighbor' is not set, and the CVs reference "
                        "more than one distinct water via nearest_water_to "
                        "-- set 'qmwater_neighbor' explicitly to disambiguate "
                        "which water belongs in the QM region."
                    )
                water_ref = next(iter(distinct_cv_refs.values()))
            elif distinct_cv_refs and _canon_water_ref(water_ref) not in distinct_cv_refs:
                raise ValueError(
                    f"'qmwater_neighbor' ({water_ref!r}) does not match any "
                    f"CV's nearest_water_to reference "
                    f"({list(distinct_cv_refs.values())!r}) -- the QM region "
                    "would then include a different water than the one used "
                    "by the collective variables / H10 mass patch. Set "
                    "'qmwater_neighbor' to match, or omit it to infer it "
                    "automatically from the CVs."
                )
            qmwater_exclude = list(inh_cfg.get("qmwater_exclude") or [])
            exclude_neighbor = inh_cfg.get("qmwater_exclude_neighbor")
            if exclude_neighbor:
                excl_resid = _find_nearest_water_resid(
                    exclude_neighbor, chain_map, top_atoms, top_residues,
                    substrate_chain_map, chain, rst7_coords,
                )
                qmwater_exclude.append(excl_resid + 1)
            resid = _find_nearest_water_resid(
                water_ref, chain_map, top_atoms, top_residues,
                substrate_chain_map, chain, rst7_coords,
                exclude_water=qmwater_exclude,
            )
            qmmask = qmmask.replace(_QMMASK_WATER_PLACEHOLDER, str(resid + 1))
        return qmmask, int(inh_cfg["qmcharge"])

    print("  Auto-selecting QM region from CV seed residues ...")
    qm_set, qmcharge, qmmask = _auto_select_qm(cv_indices_1based, parm_path)
    n = len(qm_set)
    print(f"  QM region: {n} atoms | estimated net charge: {qmcharge}")
    print(f"  qmmask: {qmmask}")
    return qmmask, qmcharge


def _auto_select_qm(
    cv_indices_1based: list[int],
    parm_path:         Path,
) -> tuple[set[int], int, str]:
    try:
        import parmed as pmd
    except ImportError:
        raise ImportError(
            "parmed is required for automatic QM region selection: pip install parmed"
        )

    top = pmd.load_file(str(parm_path))

    # Build adjacency list (0-based atom indices)
    adj: dict[int, set[int]] = defaultdict(set)
    for bond in top.bonds:
        i, j = bond.atom1.idx, bond.atom2.idx
        adj[i].add(j)
        adj[j].add(i)

    # Seed: all residues containing a CV atom
    seed_0        = {i - 1 for i in cv_indices_1based}
    seed_residues = {top.atoms[i].residue.idx for i in seed_0}

    # Start with all atoms of seed residues
    qm_set: set[int] = set()
    for res_idx in seed_residues:
        for atom in top.residues[res_idx].atoms:
            qm_set.add(atom.idx)

    qm_set   = _expand_to_cc_cuts(qm_set, adj, top)
    qmcharge = _estimate_charge(qm_set, top)
    qmmask   = _build_mask(qm_set)

    return qm_set, qmcharge, qmmask


def _element(atom) -> str:
    """Return uppercase element symbol for a parmed Atom."""
    try:
        el = atom.element
        if el:
            return str(el).upper()
    except AttributeError:
        pass
    # Fallback: first alphabetic character of atom name
    name = atom.name.lstrip("0123456789")
    return name[0].upper() if name else "?"


def _is_carbon(atom) -> bool:
    return _element(atom) == "C"


def _is_nitrogen(atom) -> bool:
    return _element(atom) == "N"


def _is_protein(res) -> bool:
    return res.name.upper() in _PROTEIN_RESNAMES


def _expand_to_cc_cuts(
    qm_set: set[int],
    adj:    dict[int, set[int]],
    top,
) -> set[int]:
    """
    Iteratively adjust the QM set until every bond crossing the QM/MM boundary
    is a C-C bond.

    For protein backbone peptide bonds (C-N):
      - If carbonyl C is QM and the amide N is MM: remove C and O of that
        residue from QM, making the CA-C bond the new boundary (C-C).
      - If amide N is QM and carbonyl C is MM: add C and O of the MM-side
        residue to QM, making the C-CA bond the new boundary (C-C).
    For all other non-C-C bonds: include the whole MM-side residue.

    Converges in at most len(residues) iterations.
    """
    for _ in range(len(top.residues) + 1):
        changed = False

        for ai in list(qm_set):
            for aj in adj.get(ai, set()):
                if aj in qm_set:
                    continue

                atom_i = top.atoms[ai]   # QM side
                atom_j = top.atoms[aj]   # MM side

                if _is_carbon(atom_i) and _is_carbon(atom_j):
                    continue  # valid C-C cut

                res_i = atom_i.residue
                res_j = atom_j.residue

                # Case A: carbonyl C (QM) -- amide N (MM): peptide C-N bond.
                # Fix: remove C and O of res_i from QM so cut becomes CA-C.
                if (atom_i.name == "C" and _is_carbon(atom_i)
                        and atom_j.name in ("N", "H", "HN") and _is_protein(res_i)):
                    for a in res_i.atoms:
                        if a.name in ("C", "O") and a.idx in qm_set:
                            qm_set.discard(a.idx)
                            changed = True

                # Case B: amide N (QM) -- carbonyl C (MM): peptide N-C bond.
                # Fix: add C and O of res_j to QM so cut becomes C-CA.
                elif (atom_i.name in ("N", "H", "HN")
                        and atom_j.name == "C" and _is_carbon(atom_j)
                        and _is_protein(res_j)):
                    for a in res_j.atoms:
                        if a.name in ("C", "O") and a.idx not in qm_set:
                            qm_set.add(a.idx)
                            changed = True

                # Case C: any other non-C-C bond -- include whole MM residue.
                else:
                    for a in res_j.atoms:
                        if a.idx not in qm_set:
                            qm_set.add(a.idx)
                            changed = True

        if not changed:
            break

    return qm_set


def _estimate_charge(qm_set: set[int], top) -> int:
    """Sum partial charges of QM atoms (parmed stores them in electron units) and round."""
    total = sum(top.atoms[i].charge for i in qm_set)
    return int(round(float(total)))


def _build_mask(qm_set: set[int]) -> str:
    """Build AMBER @-style atom mask from 0-based atom index set."""
    if not qm_set:
        return ""
    indices = sorted(i + 1 for i in qm_set)  # 1-based
    ranges  = []
    start   = end = indices[0]
    for idx in indices[1:]:
        if idx == end + 1:
            end = idx
        else:
            ranges.append(f"{start}" if start == end else f"{start}-{end}")
            start = end = idx
    ranges.append(f"{start}" if start == end else f"{start}-{end}")
    return "@" + ",".join(ranges)


# ---------------------------------------------------------------------------
# H10 hydrogen identification
# ---------------------------------------------------------------------------

def find_h10_atoms(
    cv_indices_1based: list[int],
    top_atoms:         list,
    top_residues:      list,
) -> dict[int, list[int]]:
    """
    Return {amber_resid_1based: [atom_0based_index, ...]} for all H atoms
    referenced by CVs.  For WAT residues both H atoms are included.
    """
    result: dict[int, list[int]] = {}

    for idx_1 in cv_indices_1based:
        atom = top_atoms[idx_1 - 1]
        if not _name_is_hydrogen(atom.name):
            continue

        res     = top_residues[atom.resid]
        res_key = atom.resid + 1  # 1-based

        if res.name == "WAT":
            h_indices = [
                a.index for a in top_atoms
                if a.resid == atom.resid and _name_is_hydrogen(a.name)
            ]
        else:
            h_indices = [atom.index]

        bucket = result.setdefault(res_key, [])
        for hi in h_indices:
            if hi not in bucket:
                bucket.append(hi)

    return result


def _name_is_hydrogen(name: str) -> bool:
    stripped = name.lstrip("0123456789")
    return bool(stripped) and stripped[0].upper() == "H"


# ---------------------------------------------------------------------------
# Protein residue count helper (used for centering mask)
# ---------------------------------------------------------------------------

def count_protein_residues(top_residues: list) -> int:
    """Return the 1-based index of the last protein residue (for cpptraj mask)."""
    last = 0
    for res in top_residues:
        if res.name.upper() in _PROTEIN_RESNAMES:
            last = res.index + 1
    return last

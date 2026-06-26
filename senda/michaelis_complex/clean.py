"""
senda.michaelis_complex.clean
Clean raw crystal-structure PDB lines for Michaelis complex preparation.

Operations (in order):
  1. Strip CONECT and END records.
  2. Strip explicitly excluded HETATM residues (default: GOL).
  3. Keep only the specified chains for ATOM records.
  4. Keep specified chains for HETATM non-water; keep HOH on any chain.
  5. Apply residue_renames (e.g. "216" → "LER").
  6. Apply protonation states from a reference structure map.
"""

from __future__ import annotations

# Amino acid protonation variant groups.
# A rename is applied only when both residue names belong to the same group.
_PROTONATION_GROUPS: tuple[frozenset[str], ...] = (
    frozenset({"HIS", "HID", "HIE", "HIP"}),
    frozenset({"CYS", "CYM", "CYX"}),
    frozenset({"ASP", "ASH"}),
    frozenset({"GLU", "GLH"}),
    frozenset({"LYS", "LYN"}),
    frozenset({"TYR", "TYM"}),
)

# Maps each variant name to its protonation group for O(1) lookup.
_GROUP_OF: dict[str, frozenset[str]] = {
    name: group
    for group in _PROTONATION_GROUPS
    for name in group
}


def _resname(line: str) -> str:
    return line[17:20].strip()


def _chain(line: str) -> str:
    return line[21]


def _resnum(line: str) -> int | None:
    try:
        return int(line[22:26])
    except ValueError:
        return None


def _rename_res(line: str, new_name: str) -> str:
    return line[:17] + f"{new_name:<3s}" + line[20:]


def build_protonation_map(
    ref_lines: list[str],
    chains:    frozenset[str],
) -> dict[tuple[str, int], str]:
    """
    Build a (chain, resnum) → resname map from ATOM records in a reference PDB.

    Only protein ATOM records are considered; HETATM (ligands, water, ions) are
    ignored, so a holo structure is safe to pass directly.
    """
    prot_map: dict[tuple[str, int], str] = {}
    for line in ref_lines:
        if not line.startswith("ATOM  "):
            continue
        ch = line[21]
        if ch not in chains:
            continue
        try:
            rn = int(line[22:26])
        except ValueError:
            continue
        prot_map[(ch, rn)] = line[17:20].strip()
    return prot_map


def clean(
    lines:           list[str],
    prot_map:        dict[tuple[str, int], str] | None = None,
    residue_renames: list[tuple[str, str]] | None = None,
    chains:          frozenset[str] = frozenset(("A", "B")),
    strip_resnames:  set[str] | None = None,
) -> list[str]:
    """
    Return cleaned PDB lines from a raw crystal structure.

    Parameters
    ----------
    lines           : raw PDB lines (with newlines)
    prot_map        : (chain, resnum) → target resname built from a reference PDB
                      via :func:`build_protonation_map`; rename applied only when
                      raw and reference names share a protonation group
    residue_renames : [(from_name, to_name), …] applied before protonation rename
    chains          : chain IDs to retain for ATOM/HETATM records
    strip_resnames  : HETATM residue names to always discard; defaults to {"GOL"}
    """
    _renames: dict[str, str] = dict(residue_renames or [])
    _strip:   set[str]       = strip_resnames if strip_resnames is not None else {"GOL"}

    out: list[str] = []

    for line in lines:
        rec = line[:6]

        if line.startswith(("CONECT", "TER", "END")):
            continue

        if rec in ("ATOM  ", "HETATM"):
            ch      = _chain(line)
            resname = _resname(line)

            if rec == "ATOM  " and ch not in chains:
                continue

            if rec == "HETATM":
                if resname in _strip:
                    continue
                if resname != "HOH" and ch not in chains:
                    continue

            # Step 5: apply residue_renames
            if resname in _renames:
                line    = _rename_res(line, _renames[resname])
                resname = _renames[resname]

            # Step 6: apply protonation from reference
            if prot_map is not None:
                rn = _resnum(line)
                if rn is not None:
                    ref_name = prot_map.get((ch, rn))
                    if ref_name is not None and ref_name != resname:
                        group = _GROUP_OF.get(resname)
                        if group is not None and ref_name in group:
                            line    = _rename_res(line, ref_name)
                            resname = ref_name

            out.append(line)
            continue

        if rec == "ANISOU":
            ch      = _chain(line)
            resname = _resname(line)

            if resname == "HOH":
                out.append(line)
                continue

            if resname in _strip or ch not in chains:
                continue

            if resname in _renames:
                line    = _rename_res(line, _renames[resname])
                resname = _renames[resname]

            if prot_map is not None:
                rn = _resnum(line)
                if rn is not None:
                    ref_name = prot_map.get((ch, rn))
                    if ref_name is not None and ref_name != resname:
                        group = _GROUP_OF.get(resname)
                        if group is not None and ref_name in group:
                            line    = _rename_res(line, ref_name)

            out.append(line)
            continue

        # Everything else (REMARK, HEADER, SEQRES, …)
        out.append(line)

    return out

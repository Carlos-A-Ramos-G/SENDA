"""
senda.michaelis_complex.clean
Clean raw crystal-structure PDB lines for Michaelis complex preparation.

Operations (in order):
  1. Strip CONECT and END records (we write our own TER later).
  2. Strip GOL residues (any chain).
  3. For ATOM: keep only chains A and B.
  4. For HETATM: keep chains A and B for non-water; keep HOH on any chain.
  5. Mirror the above chain rules for ANISOU records.
  6. Rename residue "216" → "LER" in ATOM/HETATM/ANISOU.
  7. Rename HIS to HIE/HID by residue number (per his_rename map).
"""

_DEFAULT_HIS: dict[int, str] = {
    41:  "HID",
    64:  "HIE",
    80:  "HID",
    163: "HIE",
    164: "HIE",
    172: "HIE",
    246: "HIE",
}

_PROTEIN_CHAINS = frozenset(("A", "B"))


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


def clean(lines: list[str],
          his_rename: dict[int, str] | None = None) -> list[str]:
    """
    Return cleaned PDB lines from a raw crystal structure.

    Parameters
    ----------
    lines      : raw PDB lines (with newlines)
    his_rename : resnum → new residue name; defaults to Mpro protonation states
    """
    if his_rename is None:
        his_rename = _DEFAULT_HIS

    out: list[str] = []

    for line in lines:
        rec = line[:6]

        # ── Records to always strip ───────────────────────────────────────────
        if line.startswith(("CONECT", "TER", "END")):
            continue

        # ── ATOM / HETATM ─────────────────────────────────────────────────────
        if rec in ("ATOM  ", "HETATM"):
            ch      = _chain(line)
            resname = _resname(line)

            if rec == "ATOM  " and ch not in _PROTEIN_CHAINS:
                continue

            if rec == "HETATM":
                if resname == "GOL":
                    continue
                if resname != "HOH" and ch not in _PROTEIN_CHAINS:
                    continue

            if resname == "216":
                line    = _rename_res(line, "LER")
                resname = "LER"

            if resname == "HIS":
                rn = _resnum(line)
                if rn is not None and rn in his_rename:
                    line = _rename_res(line, his_rename[rn])

            out.append(line)
            continue

        # ── ANISOU ────────────────────────────────────────────────────────────
        if rec == "ANISOU":
            ch      = _chain(line)
            resname = _resname(line)

            if resname == "HOH":
                out.append(line)
                continue

            if resname == "GOL" or ch not in _PROTEIN_CHAINS:
                continue

            if resname == "216":
                line    = _rename_res(line, "LER")
                resname = "LER"

            if resname == "HIS":
                rn = _resnum(line)
                if rn is not None and rn in his_rename:
                    line = _rename_res(line, his_rename[rn])

            out.append(line)
            continue

        # ── Everything else (REMARK, HEADER, SEQRES, …) ──────────────────────
        out.append(line)

    return out

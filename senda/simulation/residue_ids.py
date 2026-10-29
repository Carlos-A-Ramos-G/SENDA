"""
senda.simulation.residue_ids
PDB residue identity (chain ID + residue number) stored in the AMBER topology.

tleap numbers residues 1..N and drops chain IDs, but configs name residues by
PDB chain and number. The mapping is built once, by pairing the non-solvent
residues of the complex PDB with the topology in order, and checked three ways:

  count     same number of non-solvent residues
  name      same residue name pair by pair (protonation variants allowed:
            HIS/HID/HIE/HIP, CYS/CYX/CYM, ASP/ASH, GLU/GLH, LYS/LYN; such
            differences are listed in a note, since the topology decides the
            protonation)
  position  one heavy atom per residue (CA, else the first heavy atom both
            share) at the same place in the PDB and in the topology's
            coordinates, after removing the overall shift (tleap solvatebox
            moves the solute) and periodic images; tolerance POSITION_TOL

It is stored as the standard %FLAG RESIDUE_NUMBER / RESIDUE_CHAINID sections
(read by ParmEd, cpptraj, pytraj and VMD; other programs ignore them), added
to the parm7 text: no other part of the file is rewritten. Solvent and ions
get chain SOLVENT_CHAIN ("_"; readers drop blank chain fields) and keep their
topology numbering.

  senda-sim annotate                                every replica's 00_prep
  python -m senda.simulation.residue_ids PARM7 RST7 PDB [PARM7 ...]
"""
from __future__ import annotations

import math
from pathlib import Path

from senda.analysis.distances import _SOLVENT_RES

POSITION_TOL  = 0.5   # Angstrom
SOLVENT_CHAIN = "_"

_VARIANTS = [
    {"HIS", "HID", "HIE", "HIP"}, {"CYS", "CYX", "CYM"}, {"ASP", "ASH"},
    {"GLU", "GLH"}, {"LYS", "LYN"},
]


def _same_name(a: str, b: str) -> bool:
    return a == b or any(a in v and b in v for v in _VARIANTS)


# ---------------------------------------------------------------------------
# parm7 / rst7 / PDB readers (plain text, no ParmEd needed)
# ---------------------------------------------------------------------------

def _parm7_sections(text: str) -> dict:
    """{FLAG: (format, [raw data lines])}"""
    out, flag, fmt = {}, None, None
    for line in text.splitlines():
        if line.startswith("%FLAG"):
            flag = line.split()[1]
            out[flag] = [None, []]
        elif line.startswith("%FORMAT") and flag:
            out[flag][0] = line[line.index("(") + 1:line.rindex(")")]
        elif flag and not line.startswith("%"):
            out[flag][1].append(line)
    return out


def _parse_section(fmt: str, lines: list) -> list:
    # e.g. 20a4, 10I8, 5E16.8
    i = 0
    while fmt[i].isdigit():
        i += 1
    kind, width = fmt[i].upper(), int(fmt[i + 1:].split(".")[0])
    vals = []
    for line in lines:
        for k in range(0, len(line), width):
            field = line[k:k + width]
            if kind == "A":
                vals.append(field.strip())
            elif field.strip():
                vals.append(int(field) if kind == "I" else float(field))
    return vals


def read_parm7_residues(parm7: Path) -> dict:
    """residue labels, first-atom pointers (0-based), atom names, atomic numbers, stored IDs (or None)"""
    sec = _parm7_sections(Path(parm7).read_text())
    get = lambda f: _parse_section(*sec[f]) if f in sec else None
    natom = get("POINTERS")[0]
    labels = get("RESIDUE_LABEL")
    ptrs = [p - 1 for p in get("RESIDUE_POINTER")] + [natom]
    labels = labels[:len(ptrs) - 1]
    numbers, chains = get("RESIDUE_NUMBER"), get("RESIDUE_CHAINID")
    ids = None
    if numbers is not None and chains is not None:
        chains = (chains + [""] * len(labels))[:len(labels)]
        ids = list(zip(chains, numbers[:len(labels)]))
    return {"labels": labels, "ptrs": ptrs, "names": get("ATOM_NAME")[:natom],
            "anum": get("ATOMIC_NUMBER"), "natom": natom, "ids": ids}


def _read_rst7(rst7: Path, natom: int) -> tuple:
    lines = Path(rst7).read_text().splitlines()
    vals = []
    for line in lines[2:]:
        vals += [float(line[k:k + 12]) for k in range(0, len(line.rstrip()), 12)]
    xyz = [vals[3 * i:3 * i + 3] for i in range(natom)]
    box = None
    rest = vals[3 * natom:]
    if len(rest) in (6, 3 * natom + 6):
        box = rest[-6:-3]
    return xyz, box


def _parse_pdb(pdb: Path) -> list:
    """non-solvent residues in file order: [(chain, number, name, {atom: xyz})]"""
    res, key = [], None
    for line in Path(pdb).read_text().splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        name = line[17:20].strip()
        if name in _SOLVENT_RES:          # before the number: hybrid-36 waters past 9999
            continue
        k = (line[21].strip() or "A", int(line[22:26]), line[26], name)
        if k != key:
            res.append((k[0], k[1], name, {}))
            key = k
        res[-1][3].setdefault(line[12:16].strip(),
                              (float(line[30:38]), float(line[38:46]), float(line[46:54])))
    return res


# ---------------------------------------------------------------------------
# Matching with checks
# ---------------------------------------------------------------------------

def match_residues(parm7: Path, pdb: Path, rst7: Path | None = None) -> list:
    """
    [(chain, number)] for every topology residue (solvent: (SOLVENT_CHAIN, topology number)).
    Raises ValueError naming the first residue that fails the count, name or
    position check. Position is checked only when *rst7* is given.
    """
    top = read_parm7_residues(parm7)
    pdb_res = _parse_pdb(pdb)
    names = {r[2] for r in pdb_res}
    known = lambda lab: any(_same_name(lab, n) for n in names)
    solute = [i for i, lab in enumerate(top["labels"]) if lab not in _SOLVENT_RES and known(lab)]

    if len(solute) != len(pdb_res):
        raise ValueError(
            f"residue count mismatch: {Path(pdb).name} has {len(pdb_res)} non-solvent residues, "
            f"{Path(parm7).name} has {len(solute)} matching ones -- wrong PDB, or residues "
            "added/removed after it was written")
    protonation = []
    for (ch, num, nm, _), i in zip(pdb_res, solute):
        if not _same_name(nm, top["labels"][i]):
            raise ValueError(
                f"residue name mismatch at topology residue {i + 1}: PDB {ch}/{nm} {num} vs "
                f"topology {top['labels'][i]} -- wrong PDB or residue order")
        if nm != top["labels"][i]:
            protonation.append(f"{ch}/{nm} {num} -> {top['labels'][i]}")
    if protonation:
        shown = ", ".join(protonation[:10]) + (f", ... (+{len(protonation) - 10})" if len(protonation) > 10 else "")
        print(f"  Note: {len(protonation)} residue(s) differ from {Path(pdb).name} only in protonation "
              f"(PDB -> topology): {shown}")

    if rst7 is not None and Path(rst7).exists():
        xyz, box = _read_rst7(rst7, top["natom"])
        pairs = []
        for (ch, num, nm, atoms), i in zip(pdb_res, solute):
            tnames = top["names"][top["ptrs"][i]:top["ptrs"][i + 1]]
            heavy = [a for a in tnames if not a.startswith("H")]
            pick = "CA" if ("CA" in atoms and "CA" in tnames and nm not in ("CA",)) else \
                   next((a for a in heavy if a in atoms), None)
            if pick is not None:
                pairs.append((ch, num, nm, i, atoms[pick], xyz[top["ptrs"][i] + tnames.index(pick)]))
        if pairs:
            shift = [sorted(t[k] - p[k] for *_, p, t in pairs)[len(pairs) // 2] for k in range(3)]
            for ch, num, nm, i, p, t in pairs:
                d = [t[k] - p[k] - shift[k] for k in range(3)]
                if box:
                    d = [d[k] - box[k] * round(d[k] / box[k]) for k in range(3)]
                dev = math.sqrt(sum(x * x for x in d))
                if dev > POSITION_TOL:
                    raise ValueError(
                        f"residue position mismatch at topology residue {i + 1}: PDB {ch}/{nm} {num} "
                        f"is {dev:.1f} A from its topology counterpart {top['labels'][i]} after "
                        "removing the overall shift -- wrong PDB, or same-named residues in a "
                        "different order")

    ids = [(SOLVENT_CHAIN, i + 1) for i in range(len(top["labels"]))]
    for (ch, num, _, _), i in zip(pdb_res, solute):
        ids[i] = (ch, num)
    return ids


# ---------------------------------------------------------------------------
# Store / read
# ---------------------------------------------------------------------------

def _format(values: list, per_line: int, width: int, kind: str) -> str:
    out = []
    for k in range(0, len(values), per_line):
        chunk = values[k:k + per_line]
        out.append("".join(f"{v:<{width}}" if kind == "a" else f"{v:>{width}d}" for v in chunk))
    return "\n".join(out) if out else ""


def write_residue_ids(parm7: Path, ids: list) -> bool:
    """Add the two flags to *parm7*; returns False if identical ones are already there."""
    top = read_parm7_residues(parm7)
    if top["ids"] is not None:
        if top["ids"] == ids:
            return False
        raise ValueError(f"{parm7} already has different residue IDs; not overwritten")
    text = Path(parm7).read_text()
    if not text.endswith("\n"):
        text += "\n"
    text += ("%FLAG RESIDUE_NUMBER\n%FORMAT(10I8)\n" + _format([n for _, n in ids], 10, 8, "I") + "\n"
             "%FLAG RESIDUE_CHAINID\n%FORMAT(20a4)\n" + _format([c for c, _ in ids], 20, 4, "a") + "\n")
    Path(parm7).write_text(text)
    return True


def annotate(parm7s: list, rst7: Path, pdb: Path) -> list:
    """Match against the first parm7 (with checks) and store the IDs in every parm7 given."""
    ids = match_residues(parm7s[0], pdb, rst7)
    written = []
    for p in parm7s:
        if read_parm7_residues(p)["labels"] != read_parm7_residues(parm7s[0])["labels"]:
            raise ValueError(f"{p} and {parm7s[0]} have different residues")
        if write_residue_ids(p, ids):
            written.append(Path(p))
    return written


def annotate_prep(prep_dir: Path, pdb: Path) -> list:
    """Annotate structure.parm7 (and structure_HMR.parm7 if present) in a 00_prep directory."""
    prep_dir = Path(prep_dir)
    parm7s = [prep_dir / n for n in ("structure.parm7", "structure_HMR.parm7") if (prep_dir / n).exists()]
    if not parm7s:
        raise FileNotFoundError(f"no structure.parm7 in {prep_dir}")
    return annotate(parm7s, prep_dir / "structure.rst7", pdb)


def chain_map_from_ids(ids: list) -> dict:
    """(chain, pdb_resnum) -> 1-based topology residue; the first residue wins a shared key."""
    out = {}
    for i, (ch, num) in enumerate(ids):
        if ch and ch != SOLVENT_CHAIN:
            out.setdefault((ch, num), i + 1)
    return out


def chains_from_ids(ids: list, labels: list) -> list:
    """Protein chains in topology order."""
    from senda.config import AMINO_ACIDS

    out = []
    for (ch, _), lab in zip(ids, labels):
        if ch and ch != SOLVENT_CHAIN and lab in AMINO_ACIDS and ch not in out:
            out.append(ch)
    return out


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 4:
        sys.exit(__doc__)
    parm7s, rst7, pdb = [sys.argv[1]] + sys.argv[4:], sys.argv[2], sys.argv[3]
    done = annotate([Path(p) for p in parm7s], Path(rst7), Path(pdb))
    print("residue IDs written to: " + (", ".join(map(str, done)) if done else "none (already present)"))

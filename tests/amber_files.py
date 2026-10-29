"""Tiny AMBER/PDB writers for tests: a consistent parm7 + rst7 + PDB without AmberTools."""

from pathlib import Path


def write_system(stem: Path, residues, ncrd=None, pdb_chains=None, pdb_order=None, shift=(0.0, 0.0, 0.0),
                 box=(50.0, 50.0, 50.0)):
    """
    residues: [(resname, [(atom, element_Z, (x, y, z)), ...]), ...] in topology order.
    Writes <stem>.prmtop/.inpcrd/.pdb. The PDB lists the residues in *pdb_order* (default: same
    order), with chain IDs *pdb_chains* (default 'A') and numbers 1..N, shifted by *shift*.
    """
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    atoms = [(a, z, xyz) for _, at in residues for a, z, xyz in at]
    natom, nres = len(atoms), len(residues)
    ptr, k = [], 1
    for _, at in residues:
        ptr.append(k)
        k += len(at)

    def block(vals, per, width, kind):
        out = []
        for i in range(0, len(vals), per):
            out.append("".join(f"{v:<{width}}" if kind == "a" else f"{v:>{width}d}" for v in vals[i:i + per]))
        return "\n".join(out) + "\n"

    pointers = [natom] + [0] * 10 + [nres] + [0] * 19
    (stem.parent / f"{stem.name}.prmtop").write_text(
        "%VERSION\n%FLAG TITLE\n%FORMAT(20a4)\n\n"
        "%FLAG POINTERS\n%FORMAT(10I8)\n" + block(pointers, 10, 8, "I")
        + "%FLAG ATOM_NAME\n%FORMAT(20a4)\n" + block([a for a, _, _ in atoms], 20, 4, "a")
        + "%FLAG RESIDUE_LABEL\n%FORMAT(20a4)\n" + block([r for r, _ in residues], 20, 4, "a")
        + "%FLAG RESIDUE_POINTER\n%FORMAT(10I8)\n" + block(ptr, 10, 8, "I")
        + "%FLAG ATOMIC_NUMBER\n%FORMAT(10I8)\n" + block([z for _, z, _ in atoms], 10, 8, "I"))
    xyz = [c for _, _, p in atoms for c in p]
    lines = [f"default_name\n{(ncrd or natom):6d}\n"]
    vals = xyz + list(box) + [90.0, 90.0, 90.0]
    for i in range(0, len(vals), 6):
        lines.append("".join(f"{v:12.7f}" for v in vals[i:i + 6]) + "\n")
    (stem.parent / f"{stem.name}.inpcrd").write_text("".join(lines))
    order = pdb_order or list(range(nres))
    chains = pdb_chains or ["A"] * nres
    pdb, n = [], 1
    for j, i in enumerate(order):
        rname, at = residues[i]
        rec = "ATOM  " if rname not in ("LIG", "WAT") else "HETATM"
        for a, _, p in at:
            q = [p[c] - shift[c] for c in range(3)]
            pdb.append(f"{rec}{n:5d} {a:<4} {rname:>3} {chains[j]}{j + 1:4d}    {q[0]:8.3f}{q[1]:8.3f}{q[2]:8.3f}  1.00  0.00\n")
            n += 1
    (stem.parent / f"{stem.name}.pdb").write_text("".join(pdb) + "END\n")
    return stem


def ala(x):
    return ("ALA", [("N", 7, (x, 0.0, 0.0)), ("CA", 6, (x + 1.5, 0.0, 0.0)), ("C", 6, (x + 2.5, 1.0, 0.0))])

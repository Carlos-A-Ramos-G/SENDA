"""
ligand_parameters.resp
Generate RESP charge-fitting input files (resp.in and resp.qin) for ligands.

For ligands there are no backbone atoms to fix — all atoms are free to vary
in the RESP fit.  Symmetry-equivalent atoms (e.g. three methyl hydrogens,
two equivalent aromatic carbons) are detected and constrained to be equal.

Symmetry detection uses RDKit canonical ranking (full graph equivalence)
with a geometry-only fallback (H atoms on the same heavy atom).

Reference: Bayly et al., J. Phys. Chem. 97, 10269 (1993).
"""

import numpy as np
from pathlib import Path
from .gaussian import _elem


_SOLVENT = {'HOH', 'WAT', 'SOL', 'TIP', 'T3P', 'Na+', 'Cl-', 'NA', 'CL'}

ATOMIC_NUMBERS = {
    'H': 1,  'C': 6,  'N': 7,  'O': 8,
    'F': 9,  'P': 15, 'S': 16, 'Cl': 17, 'Br': 35,
    'Mg': 12, 'Zn': 30, 'Fe': 26, 'Mn': 25,
    'Ni': 28, 'Cu': 29, 'Co': 27, 'Se': 34,
}


def _parse_pdb(path):
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if line[17:20].strip() in _SOLVENT:
                continue
            try:
                atoms.append({
                    'name': line[12:16].strip(),
                    'x':    float(line[30:38]),
                    'y':    float(line[38:46]),
                    'z':    float(line[46:54]),
                })
            except ValueError:
                continue
    return atoms


def _pos(a):
    return np.array([a['x'], a['y'], a['z']])


def _find_equiv_rdkit(atoms, all_indices):
    """
    Symmetry-equivalent atoms via RDKit canonical ranking.

    Builds the molecular graph from covalent-radii distances, runs the Morgan
    algorithm with breakTies=False, and groups atoms that share the same rank.
    Returns {global_1based: ref_global_1based} for constrained atoms.
    """
    from rdkit import Chem

    mol = Chem.RWMol()
    for a in atoms:
        mol.AddAtom(Chem.Atom(_elem(a['name'])))

    conf = Chem.Conformer(len(atoms))
    for i, a in enumerate(atoms):
        p = _pos(a)
        conf.SetAtomPosition(i, (float(p[0]), float(p[1]), float(p[2])))
    mol.AddConformer(conf, assignId=True)

    pt = Chem.GetPeriodicTable()
    n  = mol.GetNumAtoms()
    for i in range(n):
        ri = pt.GetRcovalent(_elem(atoms[i]['name']))
        for j in range(i + 1, n):
            rj = pt.GetRcovalent(_elem(atoms[j]['name']))
            if np.linalg.norm(_pos(atoms[i]) - _pos(atoms[j])) < 1.3 * (ri + rj):
                mol.AddBond(i, j, Chem.BondType.SINGLE)

    # Skip valence check — unusual charge states don't abort sanitization
    Chem.SanitizeMol(
        mol,
        Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES,
    )

    ranks = list(Chem.CanonicalRankAtoms(mol, breakTies=False))

    rank_groups = {}
    for i in all_indices:
        rank_groups.setdefault(ranks[i], []).append(i)

    equiv = {}
    for idxs in rank_groups.values():
        if len(idxs) < 2:
            continue
        idxs = sorted(idxs)
        ref  = idxs[0] + 1      # 1-based
        for later in idxs[1:]:
            equiv[later + 1] = ref
    return equiv


def _find_equiv_geom(atoms, all_indices):
    """
    Geometry-only fallback: constrain H atoms bonded to the same heavy atom.
    Returns {global_1based: ref_global_1based}.
    """
    h_to_heavy = {}
    for li in all_indices:
        a = atoms[li]
        if _elem(a['name']) != 'H':
            continue
        for lj in all_indices:
            if li == lj or _elem(atoms[lj]['name']) not in ('C', 'N', 'O', 'S', 'P'):
                continue
            if np.linalg.norm(_pos(a) - _pos(atoms[lj])) < 1.35:
                h_to_heavy[li] = lj
                break

    heavy_to_hs = {}
    for h_idx, heavy_idx in h_to_heavy.items():
        heavy_to_hs.setdefault(heavy_idx, []).append(h_idx)

    equiv = {}
    for _, hs in heavy_to_hs.items():
        hs.sort()
        ref_global = hs[0] + 1
        for later in hs[1:]:
            equiv[later + 1] = ref_global
    return equiv


def _find_equiv(atoms):
    """
    Dispatcher: RDKit (full graph symmetry, preferred) or geometry-only fallback.
    Returns {global_1based: ref_global_1based} for symmetry-constrained atoms.
    """
    all_indices = list(range(len(atoms)))
    try:
        return _find_equiv_rdkit(atoms, all_indices)
    except ImportError:
        print("  Warning: RDKit not found -- falling back to geometry-only equivalence "
              "detection (H atoms on the same heavy atom only). Install RDKit for "
              "full symmetry detection.")
        return _find_equiv_geom(atoms, all_indices)
    except Exception as e:
        print(f"  Warning: RDKit equivalence detection failed ({e}) -- "
              "falling back to geometry-only detection.")
        return _find_equiv_geom(atoms, all_indices)


def write_resp_in(pdb_path, charge, resname, output):
    """
    Write the RESP control file (resp.in) for a ligand.

    All atoms are free to vary (no backbone constraints); symmetry-equivalent
    atoms are constrained to be equal.

    Parameters
    ----------
    pdb_path : str | Path -- ligand PDB file
    charge   : int        -- net molecular charge
    resname  : str        -- ligand residue name (used as title in the file)
    output   : str | Path -- output path (e.g. 'resp.in')
    """
    atoms = _parse_pdb(pdb_path)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records in {pdb_path}")

    equiv = _find_equiv(atoms)

    table = []
    for global_1, a in enumerate(atoms, start=1):
        elem = _elem(a['name'])
        anum = ATOMIC_NUMBERS.get(elem)
        if anum is None:
            raise ValueError(
                f"Unknown element '{elem}' (atom '{a['name']}') in {pdb_path}. "
                "Add it to ATOMIC_NUMBERS in resp.py."
            )
        constraint = equiv.get(global_1, 0)   # 0 = free; N = constrained to atom N
        table.append((anum, constraint))

    natoms = len(table)
    lines = [
        "ligand-resp run #1",
        " &cntrl",
        " nmol=1,",
        " ihfree=1,",
        " qwt=0.0005,",
        " iqopt=2,",
        " /",
        "    1.00000",
        resname,
        f"{charge:5d}{natoms:5d}",
    ]
    for anum, constraint in table:
        lines.append(f"{anum:5d}{constraint:5d}")
    lines += ["", ""]   # two trailing blank lines required by resp

    Path(output).write_text('\n'.join(lines) + '\n')
    n_equiv = sum(1 for v in equiv.values())
    print(f"  resp.in  : {output}  ({natoms} atoms, {n_equiv} symmetry-constrained)")


def write_resp_qin(pdb_path, output):
    """
    Write the initial charges file (resp.qin) for a ligand.

    All starting charges are 0.0 because all atoms are free (no pre-set backbone
    charges, unlike the amino acid protocol in CRAFT).

    Parameters
    ----------
    pdb_path : str | Path -- ligand PDB file
    output   : str | Path -- output path (e.g. 'resp.qin')
    """
    atoms   = _parse_pdb(pdb_path)
    natoms  = len(atoms)
    charges = [0.0] * natoms

    lines = []
    for i in range(0, len(charges), 8):
        lines.append(''.join(f"{q:10.6f}" for q in charges[i:i + 8]))

    Path(output).write_text('\n'.join(lines) + '\n')
    print(f"  resp.qin : {output}  ({natoms} atoms, all free)")

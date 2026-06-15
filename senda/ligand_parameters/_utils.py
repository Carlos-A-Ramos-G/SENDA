"""
ligand_parameters._utils
Shared utilities: PDB residue name reading and charge detection.
"""

from pathlib import Path


_SOLVENT_RESNAMES = {"HOH", "WAT", "SOL", "TIP", "T3P", "Na+", "Cl-", "NA", "CL"}


def get_resname(pdb_path) -> str:
    """Read the residue name from the first non-solvent HETATM or ATOM record."""
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith(("HETATM", "ATOM  ")):
                resname = line[17:20].strip()
                if resname and resname not in _SOLVENT_RESNAMES:
                    return resname
    raise RuntimeError(
        f"No valid residue name found in {pdb_path}. "
        "Ensure the file contains HETATM or ATOM records with a 3-letter residue name."
    )


def infer_net_charge(pdb_path) -> int:
    """
    Infer net formal charge from a ligand PDB using RDKit.

    Tries candidate total charges until DetermineBonds produces a
    self-consistent assignment.
    """
    try:
        from rdkit import Chem
        from rdkit.Chem import rdDetermineBonds
    except ImportError:
        raise RuntimeError(
            "RDKit is required for charge detection. "
            "Install with: conda install -c conda-forge rdkit"
        )

    mol = Chem.MolFromPDBFile(str(pdb_path), removeHs=False, sanitize=False)
    if mol is None:
        raise RuntimeError(f"RDKit could not parse {pdb_path}")

    for charge in [0, -1, 1, -2, 2, -3, 3]:
        try:
            mol_try = Chem.RWMol(Chem.Mol(mol))
            rdDetermineBonds.DetermineBonds(mol_try, charge=charge)
            Chem.SanitizeMol(mol_try)
            if Chem.GetFormalCharge(mol_try) == charge:
                return charge
        except Exception:
            continue

    raise RuntimeError(
        f"Could not determine net charge for {Path(pdb_path).name}. "
        "Check that the PDB contains all hydrogens and correct geometry, "
        "or set net_charge to an explicit integer in config.yaml."
    )


def resolve_charge(pdb_path, lcfg: dict) -> int:
    """Resolve net charge: auto-detect via RDKit or use explicit value from config."""
    raw = lcfg.get("net_charge", "auto")
    if str(raw).lower() == "auto":
        return infer_net_charge(pdb_path)
    return int(raw)

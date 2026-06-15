"""
ligand_parameters.gaussian
Generate Gaussian input files (.com) from a ligand PDB or optimised log.

Two flavours:
  write_com    – geometry optimisation (B3LYP/6-31G*, reads ligand PDB directly)
  write_hf_com – HF/6-31G(d) single-point ESP (from optimised geometry)
"""

from pathlib import Path


NPROC_DEFAULT    = 16
MEM_DEFAULT      = "512MB"
ROUTE_DEFAULT    = "#P b3lyp/6-31g* opt"
HF_ROUTE_DEFAULT = "#p hf/6-31g(d) SCF=Tight Pop=MK IOp(6/33=2)"

# Atomic number → element symbol (used when parsing Gaussian opt logs)
ATOMIC_SYMBOLS = {
    1: 'H',   6: 'C',   7: 'N',   8: 'O',
    9: 'F',  15: 'P',  16: 'S',  17: 'Cl',  35: 'Br',
   12: 'Mg', 26: 'Fe', 25: 'Mn', 28: 'Ni',
   29: 'Cu', 27: 'Co', 30: 'Zn', 34: 'Se',
}

# Two-letter elements that must be detected before falling back to first char
_ELEM2 = {
    'CL': 'Cl', 'BR': 'Br', 'MG': 'Mg', 'ZN': 'Zn',
    'FE': 'Fe', 'MN': 'Mn', 'NI': 'Ni', 'CU': 'Cu',
    'CO': 'Co', 'SE': 'Se',
}

_SOLVENT = {'HOH', 'WAT', 'SOL', 'TIP', 'T3P', 'Na+', 'Cl-', 'NA', 'CL'}


def _elem(atom_name: str) -> str:
    """Extract element symbol from a PDB atom name, handling 2-letter elements."""
    name = atom_name.strip().lstrip('0123456789').upper()
    if len(name) >= 2 and name[:2] in _ELEM2:
        return _ELEM2[name[:2]]
    return name[0] if name else atom_name.strip()[0]


def _parse_pdb(path):
    """Return list of {name, x, y, z} dicts for each non-solvent ATOM/HETATM."""
    atoms = []
    with open(path) as f:
        for line in f:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            if line[17:20].strip() in _SOLVENT:
                continue
            atoms.append({
                'name': line[12:16].strip(),
                'x':    float(line[30:38]),
                'y':    float(line[38:46]),
                'z':    float(line[46:54]),
            })
    return atoms


def _write_gjf(com_path, header_lines, title, charge, mult, atoms_xyz):
    """Shared Gaussian input writer: header + title + charge/mult + coords."""
    blocks = list(header_lines) + ["", title, "", f"{charge} {mult}"]
    for elem, x, y, z in atoms_xyz:
        blocks.append(f" {elem:<2s}  {x:16.8f}{y:16.8f}{z:16.8f}")
    blocks.append("")   # Gaussian requires trailing blank line
    Path(com_path).write_text('\n'.join(blocks) + '\n')


def write_com(pdb_path, com_path, charge, mult,
              nproc=NPROC_DEFAULT, mem=MEM_DEFAULT, route=ROUTE_DEFAULT):
    """
    Write a geometry-optimisation .com from a ligand PDB.

    No capping is required for ligands — the PDB is used as-is.
    """
    atoms = _parse_pdb(pdb_path)
    if not atoms:
        raise ValueError(f"No ATOM/HETATM records in {pdb_path}")

    base  = Path(pdb_path).stem
    title = f"{base}  b3lyp/6-31g* geometry optimisation"

    header    = [f"%nprocshared={nproc}", f"%mem={mem}", route]
    atoms_xyz = [(_elem(a['name']), a['x'], a['y'], a['z']) for a in atoms]
    _write_gjf(com_path, header, title, charge, mult, atoms_xyz)

    print(f"Input  : {pdb_path}  ({len(atoms)} atoms)")
    print(f"Output : {com_path}")
    print(f"  Charge / mult : {charge} / {mult}")
    print(f"  nproc / mem   : {nproc} / {mem}")


def parse_opt_log(log_path):
    """
    Extract the final optimised geometry from a Gaussian optimisation log.

    Finds the LAST 'Standard orientation:' block and returns a list of
    (element_symbol, x, y, z) tuples in file order.
    """
    lines = Path(log_path).read_text().splitlines()

    last_block = None
    for i, line in enumerate(lines):
        if 'Standard orientation:' in line:
            last_block = i

    if last_block is None:
        raise ValueError(f"No 'Standard orientation:' block in {log_path}")

    atoms = []
    i = last_block + 5          # skip title + 3 header lines + separator
    while i < len(lines):
        line = lines[i]
        if '---' in line:
            break
        parts = line.split()
        if len(parts) == 6:
            atomic_num = int(parts[1])
            x, y, z = float(parts[3]), float(parts[4]), float(parts[5])
            elem = ATOMIC_SYMBOLS.get(atomic_num)
            if elem is None:
                raise ValueError(
                    f"Unknown atomic number {atomic_num} in {log_path}. "
                    "Add it to ATOMIC_SYMBOLS in gaussian.py."
                )
            atoms.append((elem, x, y, z))
        i += 1

    if not atoms:
        raise ValueError(
            f"Could not parse atoms from Standard orientation block in {log_path}")
    return atoms


def write_hf_com(atoms_xyz, com_path, base, charge, mult,
                 nproc=NPROC_DEFAULT, mem=MEM_DEFAULT, route=HF_ROUTE_DEFAULT):
    """
    Write an HF/6-31G(d) single-point .com for ESP/RESP from an atom list.

    atoms_xyz : list of (element, x, y, z) — typically from parse_opt_log()
    base      : ligand name stem (e.g. 'LIG'), used for the file title
    """
    header = [f"%nprocshared={nproc}", f"%mem={mem}", route]
    _write_gjf(com_path, header, f"{base}", charge, mult, atoms_xyz)

    print(f"Output : {com_path}")
    print(f"  Charge / mult : {charge} / {mult}")
    print(f"  nproc / mem   : {nproc} / {mem}")
    print(f"  Atoms         : {len(atoms_xyz)}")

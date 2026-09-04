"""
senda.qmmm.common.h10

Generate the H10 topology: set mass of QM-involved hydrogen atoms to 10 amu.
For WAT residues, both H atoms are patched even if only one is in the CVs.
"""
from pathlib import Path


def generate_h10_topology(
    parm_path:  Path,
    h10_atoms:  dict,
    out_path:   Path,
) -> None:
    """
    Write a modified parm7 file with H masses set to 10.

    h10_atoms: {amber_resid_1based: [atom_0based_index, ...]}
               as returned by atoms.find_h10_atoms().
    """
    try:
        import parmed as pmd
    except ImportError:
        raise ImportError("parmed is required for H10 topology generation: pip install parmed")

    if not h10_atoms:
        raise ValueError(
            "No hydrogen atoms found in the collective variables. "
            "Cannot generate H10 topology -- verify CV atom names."
        )

    top = pmd.load_file(str(parm_path))

    patched = []
    for res_key, atom_indices in h10_atoms.items():
        for ai in atom_indices:
            top.atoms[ai].mass = 10.0
            patched.append(f"{top.atoms[ai].name} (resid {res_key})")

    top.save(str(out_path), overwrite=True)
    print(f"  H10 topology: {out_path.name}  [{len(patched)} H atoms patched]")
    for label in patched:
        print(f"    {label}")

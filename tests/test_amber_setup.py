"""Regression tests for the AMBER-path script generation (no AMBER/OpenMM needed)."""

import re
import subprocess
import sys

from senda.simulation.templates import make_restrainer
from senda.workflow.run import _param_script

# --- restrainer.py ----------------------------------------------------------

def _fake_prmtop(residues):
    """Minimal parm7 with ATOM_NAME / RESIDUE_LABEL / RESIDUE_POINTER sections."""
    names, labels, ptrs = [], [], []
    for resname, atoms in residues:
        ptrs.append(len(names) + 1)
        labels.append(resname)
        names.extend(atoms)

    def section(flag, fmt, width, items):
        per = 80 // width
        body = "\n".join("".join(f"{x:<{width}}" if isinstance(x, str) else f"{x:>{width}}"
                                 for x in items[i:i + per]) for i in range(0, len(items), per))
        return f"%FLAG {flag}\n%FORMAT({fmt})\n{body}\n"

    return (section("ATOM_NAME", "20a4", 4, names)
            + section("RESIDUE_LABEL", "20a4", 4, labels)
            + section("RESIDUE_POINTER", "10I8", 8, ptrs))


LIG = ["C1", "C2", "C3", "C4"]
TOPOLOGY = [("ALA", ["N", "CA"]), ("LER", LIG), ("ALA", ["N", "CA"]), ("LER", LIG)]
# atom numbers (1-based): ALA 1-2, LER copy0 3-6, ALA 7-8, LER copy1 9-12

WELL = {"r1": -10.0, "r2": -5.0, "r3": 15.0, "r4": 20.0, "rk2": 500, "rk3": 500}


def _run_restrainer(tmp_path, restraints):
    (tmp_path / "00_prep").mkdir()
    (tmp_path / "00_prep" / "structure.parm7").write_text(_fake_prmtop(TOPOLOGY))
    (tmp_path / "restrainer.py").write_text(make_restrainer(restraints))
    subprocess.run([sys.executable, "restrainer.py"], cwd=tmp_path, check=True,
                   capture_output=True)
    return [list(map(int, m.group(1).split(",")[:4]))
            for m in re.finditer(r"iat=([\d,\s]+),\n", (tmp_path / "restr").read_text())]


def _dihedral(**extra):
    return {"atoms": [{"residue": "LER", "name": n, **extra} for n in LIG], **WELL}


def test_same_residue_without_index_applies_to_every_copy(tmp_path):
    assert _run_restrainer(tmp_path, [_dihedral()]) == [[3, 4, 5, 6], [9, 10, 11, 12]]


def test_same_residue_with_index_selects_one_copy(tmp_path):
    # chain: B resolves to index 1 at setup; previously ignored -> both copies
    assert _run_restrainer(tmp_path, [_dihedral(index=1)]) == [[9, 10, 11, 12]]


def test_index_on_one_atom_is_inherited_by_same_residue(tmp_path):
    r = _dihedral()
    r["atoms"][0]["index"] = 1
    assert _run_restrainer(tmp_path, [r]) == [[9, 10, 11, 12]]


def test_same_resname_across_copies(tmp_path):
    r = {"atoms": [{"residue": "LER", "name": "C1", "index": 0},
                   {"residue": "LER", "name": "C1", "index": 1}], **WELL}
    assert _run_restrainer(tmp_path, [r]) == [[3, 9, 0, 0]]


def test_cross_residue_defaults_to_first_copy(tmp_path):
    r = {"atoms": [{"residue": "ALA", "name": "CA"}, {"residue": "LER", "name": "C1", "index": 1}],
         **WELL}
    assert _run_restrainer(tmp_path, [r]) == [[2, 9, 0, 0]]


# --- senda-slurm inserts amber_module verbatim ------------------------------

def test_param_script_inserts_amber_module_verbatim():
    slurm = {"amber_module": "module load amber/24", "cpu": {"time": "1:00:00"}}
    script = _param_script(slurm, "config.yaml", "logs")
    assert "\nmodule load amber/24\n" in script
    assert "module load module load" not in script

"""Regression tests for the AMBER-path script generation (no AMBER/OpenMM needed)."""

import re
import subprocess
import sys

from senda.simulation.setup import _write_input_files, _write_run_scripts
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


# --- equilibration cycle count follows equil_schedule -----------------------

SIM = {
    "temperature": 300.0,
    "min": {"maxcyc": 100, "ncyc": 10, "cut": 10.0, "max_cycles_cap": 2,
            "convergence_threshold": 1e-3},
    "heat": {"ps": 10},
    "equil": {"npt_ns": 0.01, "nvt_ns": 0.01},
    "production": {"ns_per_chunk": 1, "total_chunks": 1, "chunks_per_job": 1},
    "output": {"ntpr": 100, "ntwx": 100, "ntwr": 100},
}


def _setup(tmp_path, schedule, mode="local"):
    for d in ("01_min", "02_heat", "03_equil", "04_NVT"):
        (tmp_path / d).mkdir(exist_ok=True)
    sim = {**SIM, "restraints": {"positional": {"equil_schedule": schedule}}}
    topology = _write_input_files(tmp_path, "APO", sim)
    _write_run_scripts(tmp_path, "APO", "WT", 1, sim, {"amber_module": "module load amber",
                       "gpu": {}}, topology, mode)
    return sim


def test_equil_cycles_follow_schedule(tmp_path):
    _setup(tmp_path, [10.0, 5.0, 2.0])
    ins = sorted(p.name for p in (tmp_path / "03_equil").glob("equil_*.in"))
    assert ins == ["equil_1.in", "equil_2.in", "equil_3.in", "equil_4.in"]
    nvt = (tmp_path / "03_equil" / "equil_4.in").read_text()
    assert "cycle 4 (unrestrained NVT)" in nvt and "ntr=1" not in nvt
    run = (tmp_path / "run_local").read_text()
    assert "for ((i=1; i<=4; i++))" in run
    assert 'cp "$DIR/03_equil/structure_equil_4.rst7" structure_NVT_0.rst7' in run


def test_default_schedule_keeps_six_cycles(tmp_path):
    _setup(tmp_path, [15.0, 12.0, 9.0, 6.0, 3.0], mode="cluster")
    run = (tmp_path / "run_gpu").read_text()
    assert "for ((i=1; i<=6; i++))" in run
    assert 'cp "$DIR/03_equil/structure_equil_6.rst7" structure_NVT_0.rst7' in run
    assert "cycle 6 (unrestrained NVT)" in (tmp_path / "03_equil" / "equil_6.in").read_text()


def test_force_resetup_removes_stale_equil_inputs(tmp_path):
    _setup(tmp_path, [15.0, 12.0, 9.0, 6.0, 3.0])
    _setup(tmp_path, [5.0])
    assert sorted(p.name for p in (tmp_path / "03_equil").glob("equil_*.in")) == \
        ["equil_1.in", "equil_2.in"]


def test_check_counts_equil_cycles_from_schedule(tmp_path):
    from senda.simulation.check import _check_replica
    (tmp_path / "00_prep").mkdir()
    (tmp_path / "00_prep" / "structure.parm7").write_text("")
    for d in ("01_min", "02_heat", "03_equil", "04_NVT"):
        (tmp_path / d).mkdir()
    done = "   Total wall time:  10    seconds\n"
    (tmp_path / "01_min" / "structure_min_1.out").write_text(done)
    (tmp_path / "02_heat" / "structure_heat.out").write_text(done)
    for i in (1, 2, 3):
        (tmp_path / "03_equil" / f"structure_equil_{i}.out").write_text(done)
    r = _check_replica(tmp_path, total_nvt_chunks=1, n_equil_amber=3)
    assert r["last_done"] == "03_equil" and r["last_done_prog"] == "3/3"



# --- senda-slurm inserts amber_module verbatim ------------------------------

def test_param_script_inserts_amber_module_verbatim():
    slurm = {"amber_module": "module load amber/24", "cpu": {"time": "1:00:00"}}
    script = _param_script(slurm, "config.yaml", "logs")
    assert "\nmodule load amber/24\n" in script
    assert "module load module load" not in script

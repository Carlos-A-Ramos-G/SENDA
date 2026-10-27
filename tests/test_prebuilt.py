"""amber_simulator.prebuilt: copy a validated topology instead of running tleap (no AMBER/OpenMM needed)."""

import pytest

from senda.simulation.setup import prebuilt_pairs, setup_all

SIM = {
    "temperature": 300.0,
    "min": {"maxcyc": 100, "ncyc": 10, "cut": 10.0, "max_cycles_cap": 2,
            "convergence_threshold": 1e-3},
    "heat": {"ps": 10},
    "equil": {"npt_ns": 0.01, "nvt_ns": 0.01},
    "production": {"ns_per_chunk": 1, "total_chunks": 1, "chunks_per_job": 1},
    "output": {"ntpr": 100, "ntwx": 100, "ntwr": 100},
}
SLURM = {"amber_module": "module load amber", "gpu": {}}


def _prebuilt(tmp_path, name, natom=3, ncrd=None):
    """<name>.prmtop/.inpcrd/.pdb with natom atoms (the inpcrd can disagree)."""
    stem = tmp_path / "builds" / name
    stem.parent.mkdir(exist_ok=True)
    pointers = "".join(f"{x:>8}" for x in [natom] + [0] * 9)
    (stem.parent / f"{name}.prmtop").write_text(
        f"%VERSION\n%FLAG TITLE\n%FORMAT(20a4)\n\n%FLAG POINTERS\n%FORMAT(10I8)\n{pointers}\n")
    (stem.parent / f"{name}.inpcrd").write_text(f"default_name\n{ncrd or natom:6d}\n")
    (stem.parent / f"{name}.pdb").write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000\nEND\n")
    return str(stem)


def _sim(tmp_path, **extra):
    return {**SIM, "prebuilt": {"L02": {"S1": _prebuilt(tmp_path, "4H4E_S1_solv")},
                                "L03": {"S1p": _prebuilt(tmp_path, "4H4E_S1p_solv")}}, **extra}


def test_prebuilt_pairs_lists_only_configured_pairs(tmp_path):
    assert prebuilt_pairs(_sim(tmp_path)) == {"L02": ["S1"], "L03": ["S1p"]}
    assert prebuilt_pairs(SIM) == {}


@pytest.mark.parametrize("engine", ["amber", "openmm"])
def test_prebuilt_copies_topology_and_skips_tleap(tmp_path, engine):
    sim = _sim(tmp_path, engine=engine)
    setup_all(["L02"], ["S1"], 1, sim, SLURM, tmp_path / "protein",
              tmp_path / "ligands_libraries", tmp_path / "simulations", mode="cluster")
    prep = tmp_path / "simulations" / "L02" / "S1" / "replica_1" / "00_prep"
    assert sorted(p.name for p in prep.iterdir() if p.suffix in (".parm7", ".rst7", ".pdb")) == \
        ["S1_L02_dimer.pdb", "structure.parm7", "structure.rst7"]
    assert not (prep / "leap_structure").exists()
    assert not any(prep.glob("*.lib")) and not any(prep.glob("*.frcmod"))
    run = (prep.parent / "run_gpu").read_text()
    assert "tleap" not in run.replace("no tleap", "")
    assert "prebuilt topology" in run and "cpptraj -i HMR.ccptraj" in run


def test_prebuilt_atom_count_mismatch_fails(tmp_path):
    sim = {**SIM, "prebuilt": {"L02": {"S1": _prebuilt(tmp_path, "bad", natom=3, ncrd=4)}}}
    with pytest.raises(ValueError, match="3 atoms"):
        setup_all(["L02"], ["S1"], 1, sim, SLURM, tmp_path / "protein",
                  tmp_path / "ligands_libraries", tmp_path / "simulations")


def test_prebuilt_missing_file_fails(tmp_path):
    sim = {**SIM, "prebuilt": {"L02": {"S1": str(tmp_path / "nowhere")}}}
    with pytest.raises(FileNotFoundError, match="L02/S1"):
        setup_all(["L02"], ["S1"], 1, sim, SLURM, tmp_path / "protein",
                  tmp_path / "ligands_libraries", tmp_path / "simulations")

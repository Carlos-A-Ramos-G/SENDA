"""amber_simulator.prebuilt: copy a validated topology instead of running tleap (no AMBER/OpenMM needed)."""

import pytest

from senda.simulation.setup import prebuilt_pairs, setup_all

from amber_files import ala, write_system

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
    """<name>.prmtop/.inpcrd/.pdb: one ALA (3 atoms); the inpcrd can claim a different atom count."""
    return str(write_system(tmp_path / "builds" / name, [ala(0.0)], ncrd=ncrd))


def _sim(tmp_path, **extra):
    return {**SIM, "prebuilt": {"INH1": {"WT": _prebuilt(tmp_path, "WT_INH1_solv")},
                                "INH2": {"MUT1": _prebuilt(tmp_path, "MUT1_INH2_solv")}}, **extra}


def test_prebuilt_pairs_lists_only_configured_pairs(tmp_path):
    assert prebuilt_pairs(_sim(tmp_path)) == {"INH1": ["WT"], "INH2": ["MUT1"]}
    assert prebuilt_pairs(SIM) == {}


@pytest.mark.parametrize("engine", ["amber", "openmm"])
def test_prebuilt_copies_topology_and_skips_tleap(tmp_path, engine):
    sim = _sim(tmp_path, engine=engine)
    setup_all(["INH1"], ["WT"], 1, sim, SLURM, tmp_path / "protein",
              tmp_path / "ligands_libraries", tmp_path / "simulations", mode="cluster")
    prep = tmp_path / "simulations" / "INH1" / "WT" / "replica_1" / "00_prep"
    assert sorted(p.name for p in prep.iterdir() if p.suffix in (".parm7", ".rst7", ".pdb")) == \
        ["WT_INH1_complex.pdb", "structure.parm7", "structure.rst7"]
    assert not (prep / "leap_structure").exists()
    assert not any(prep.glob("*.lib")) and not any(prep.glob("*.frcmod"))
    run = (prep.parent / "run_gpu").read_text()
    assert "tleap" not in run.replace("no tleap", "")
    assert "prebuilt topology" in run and "cpptraj -i HMR.ccptraj" in run
    from senda.simulation.residue_ids import read_parm7_residues
    assert read_parm7_residues(prep / "structure.parm7")["ids"] == [("A", 1)]   # IDs stored at setup


def test_prebuilt_atom_count_mismatch_fails(tmp_path):
    sim = {**SIM, "prebuilt": {"INH1": {"WT": _prebuilt(tmp_path, "bad", natom=3, ncrd=4)}}}
    with pytest.raises(ValueError, match="3 atoms"):
        setup_all(["INH1"], ["WT"], 1, sim, SLURM, tmp_path / "protein",
                  tmp_path / "ligands_libraries", tmp_path / "simulations")


def test_prebuilt_missing_file_fails(tmp_path):
    sim = {**SIM, "prebuilt": {"INH1": {"WT": str(tmp_path / "nowhere")}}}
    with pytest.raises(FileNotFoundError, match="INH1/WT"):
        setup_all(["INH1"], ["WT"], 1, sim, SLURM, tmp_path / "protein",
                  tmp_path / "ligands_libraries", tmp_path / "simulations")

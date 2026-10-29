"""Complex PDB naming ({mut}_{inh}_complex.pdb, legacy _dimer.pdb) and chain defaults for any number of chains."""

import pytest

from senda.config import (complex_pdb_name, find_complex_pdb, protein_chains, resolve_chains,
                          system_chains)


def _pdb(path, chains, ligand_chain=None):
    """One ALA per chain (blank chain allowed as ' '), an optional ligand and a hybrid-36 water."""
    lines = [f"ATOM  {i + 1:5d}  CA  ALA {ch}{i + 1:4d}       0.000   0.000   0.000  1.00  0.00           C\n"
             for i, ch in enumerate(chains)]
    if ligand_chain:
        lines.append(f"HETATM   90  C1  LIG {ligand_chain} 900       1.000   0.000   0.000  1.00  0.00           C\n")
    lines.append("ATOM     99  O   WAT  A000       2.000   0.000   0.000  1.00  0.00           O\n")
    path.write_text("".join(lines) + "END\n")
    return path


def test_complex_pdb_name():
    assert complex_pdb_name("WT", "L01") == "WT_L01_complex.pdb"


def test_find_prefers_current_name_then_legacy(tmp_path):
    a, b = tmp_path / "prep", tmp_path / "protein"
    a.mkdir(); b.mkdir()
    (b / "WT_L01_dimer.pdb").write_text("END\n")
    assert find_complex_pdb([a, b], "WT", "L01") == b / "WT_L01_dimer.pdb"
    (b / "WT_L01_complex.pdb").write_text("END\n")
    assert find_complex_pdb([a, b], "WT", "L01") == b / "WT_L01_complex.pdb"
    with pytest.raises(FileNotFoundError, match="WT_L02_complex.pdb"):
        find_complex_pdb([a, b], "WT", "L02")


@pytest.mark.parametrize("chains, expected", [
    ("A", ["A"]),                      # monomer
    (" ", ["A"]),                      # monomer without chain ID
    ("AB", ["A", "B"]),                # dimer
    ("ABCD", ["A", "B", "C", "D"]),    # tetramer
])
def test_protein_chains(tmp_path, chains, expected):
    assert protein_chains(_pdb(tmp_path / "x.pdb", chains)) == expected


def test_ligand_on_its_own_chain_is_not_a_protein_chain(tmp_path):
    assert protein_chains(_pdb(tmp_path / "x.pdb", "AB", ligand_chain="L")) == ["A", "B"]


def test_config_chains_win(tmp_path):
    pdb = _pdb(tmp_path / "x.pdb", "ABCD")
    assert resolve_chains({"michaelis_complex": {"chains": ["B"]}}, pdb) == ["B"]
    assert resolve_chains({}, pdb) == ["A", "B", "C", "D"]
    assert resolve_chains({}) == ["A"]


def test_system_chains_reads_the_replica_pdb(tmp_path):
    prep = tmp_path / "simulations" / "L01" / "WT" / "replica_1" / "00_prep"
    prep.mkdir(parents=True)
    _pdb(prep / "WT_L01_dimer.pdb", " ")              # legacy name, monomer
    assert system_chains({}, tmp_path, "L01", "WT") == ["A"]
    _pdb(prep / "WT_L01_complex.pdb", "AB")
    assert system_chains({}, tmp_path, "L01", "WT") == ["A", "B"]
    assert system_chains({}, tmp_path, "L02", "WT") == ["A"]   # nothing found

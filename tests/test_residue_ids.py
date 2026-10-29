"""Residue IDs (PDB chain + number) stored in the topology, with count/name/position checks."""

import pytest

from senda.analysis.distances import _build_chain_map
from senda.simulation.residue_ids import (SOLVENT_CHAIN, annotate, chains_from_ids, chain_map_from_ids,
                                          match_residues, read_parm7_residues)

from amber_files import ala, write_system


def lig(x):
    return ("LIG", [("C1", 6, (x, 5.0, 0.0)), ("O1", 8, (x + 1.2, 5.0, 0.0))])


WAT = ("WAT", [("O", 8, (40.0, 40.0, 40.0)), ("H1", 1, (40.9, 40.0, 40.0)), ("H2", 1, (39.7, 40.9, 40.0))])


def _files(stem):
    return stem.parent / f"{stem.name}.prmtop", stem.parent / f"{stem.name}.inpcrd", stem.parent / f"{stem.name}.pdb"


def _dimer(tmp_path, **kw):
    # topology order: protein A, ligand A, protein B, ligand B, water
    res = [ala(0.0), lig(0.0), ala(20.0), lig(20.0), WAT]
    return _files(write_system(tmp_path / "sys", res, pdb_chains=list("AABB"), pdb_order=[0, 1, 2, 3], **kw))


def _numbers(pdb):
    # renumber the PDB the way a crystal file would: protein 1 per chain, ligand 301
    text = pdb.read_text().replace(" A   2 ", " A 301 ").replace(" B   3 ", " B   1 ").replace(" B   4 ", " B 301 ")
    pdb.write_text(text)


def test_dimer_round_trip(tmp_path):
    parm7, rst7, pdb = _dimer(tmp_path)
    _numbers(pdb)
    assert annotate([parm7], rst7, pdb) == [parm7]
    top = read_parm7_residues(parm7)
    assert top["ids"] == [("A", 1), ("A", 301), ("B", 1), ("B", 301), (SOLVENT_CHAIN, 5)]
    assert chain_map_from_ids(top["ids"]) == {("A", 1): 1, ("A", 301): 2, ("B", 1): 3, ("B", 301): 4}
    assert chains_from_ids(top["ids"], top["labels"]) == ["A", "B"]
    assert annotate([parm7], rst7, pdb) == []                  # idempotent


def test_monomer_blank_chain_and_hybrid36_water(tmp_path):
    parm7, rst7, pdb = _files(write_system(tmp_path / "m", [ala(0.0), lig(0.0), WAT], pdb_chains=["A", "A"],
                                           pdb_order=[0, 1]))
    text = pdb.read_text().replace(" A   ", "     ")         # blank chain IDs
    pdb.write_text(text.replace("END\n", "HETATM 9999  O   WAT  A000      40.000  40.000  40.000  1.00  0.00\nEND\n"))
    assert match_residues(parm7, pdb, rst7)[:2] == [("A", 1), ("A", 2)]


def test_uniform_shift_is_accepted(tmp_path):
    parm7, rst7, pdb = _dimer(tmp_path, shift=(7.0, -3.0, 2.0))     # e.g. tleap solvatebox
    assert match_residues(parm7, pdb, rst7)[0] == ("A", 1)


def test_name_mismatch_is_caught(tmp_path):
    # crystal-style PDB: both proteins first, then both ligands (topology regrouped by chain)
    res = [ala(0.0), lig(0.0), ala(20.0), lig(20.0)]
    parm7, rst7, pdb = _files(write_system(tmp_path / "s", res, pdb_chains=list("ABAB"), pdb_order=[0, 2, 1, 3]))
    with pytest.raises(ValueError, match="name mismatch at topology residue 2"):
        match_residues(parm7, pdb, rst7)


def test_identical_ligands_swapped_are_caught_by_position(tmp_path):
    res = [ala(0.0), ala(20.0), lig(0.0), lig(20.0)]
    parm7, rst7, pdb = _files(write_system(tmp_path / "s", res, pdb_chains=list("ABBA"), pdb_order=[0, 1, 3, 2]))
    with pytest.raises(ValueError, match="position mismatch at topology residue 3"):
        match_residues(parm7, pdb, rst7)


def test_count_mismatch_is_caught(tmp_path):
    parm7, rst7, pdb = _dimer(tmp_path)
    pdb.write_text("".join(l for l in pdb.read_text().splitlines(True) if " B   4" not in l))
    with pytest.raises(ValueError, match="count mismatch"):
        match_residues(parm7, pdb, rst7)


def test_different_ids_are_not_overwritten(tmp_path):
    parm7, rst7, pdb = _dimer(tmp_path)
    annotate([parm7], rst7, pdb)
    _numbers(pdb)
    with pytest.raises(ValueError, match="already has different residue IDs"):
        annotate([parm7], rst7, pdb)


def test_chain_map_from_topology_needs_no_pdb(tmp_path):
    parm7, rst7, pdb = _dimer(tmp_path)
    _numbers(pdb)
    unannotated = _build_chain_map(parm7, pdb)                        # matched in memory, file untouched
    assert read_parm7_residues(parm7)["ids"] is None
    with pytest.raises(FileNotFoundError, match="senda-sim annotate"):
        _build_chain_map(parm7, None)
    annotate([parm7], rst7, pdb)
    assert _build_chain_map(parm7, None) == unannotated == {("A", 1): 1, ("A", 301): 2, ("B", 1): 3, ("B", 301): 4}


def test_protonation_differences_are_reported_not_rejected(tmp_path, capsys):
    his = ("HIP", [("N", 7, (0.0, 0.0, 0.0)), ("CA", 6, (1.5, 0.0, 0.0)), ("C", 6, (2.5, 1.0, 0.0))])
    parm7, rst7, pdb = _files(write_system(tmp_path / "h", [his, lig(0.0)], pdb_order=[0, 1]))
    pdb.write_text(pdb.read_text().replace(" HIP A   1", " HIS A  37"))
    assert match_residues(parm7, pdb, rst7)[0] == ("A", 37)
    assert "1 residue(s) differ from h.pdb only in protonation (PDB -> topology): A/HIS 37 -> HIP" in capsys.readouterr().out

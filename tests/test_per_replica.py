"""senda-analyse --per-replica: best frame of each replica against the pooled modes (no pytraj needed)."""

import numpy as np

from senda.analysis.distances import _select_frame, _select_per_replica


def _data():
    # 2 distances, 3 replicas x 4 frames; modes 3.0 / 5.0, sigma 0.5
    pooled = np.array([
        [3.9, 5.9], [3.1, 5.1], [3.6, 5.6], [4.0, 6.0],     # replica 1: best frame 2
        [3.4, 5.4], [3.3, 5.3], [3.8, 5.0], [3.0, 5.2],     # replica 2: best frame 4
        [4.5, 6.5], [4.2, 6.4], [4.4, 6.1], [4.6, 6.6],     # replica 3: nothing within 1 sigma
    ])
    lookup = [(rep, f) for rep in (1, 2, 3) for f in (1, 2, 3, 4)]
    return pooled, [3.0, 5.0], [0.5, 0.5], lookup


def test_each_replica_gets_its_own_best_frame():
    pooled, modes, stds, lookup = _data()
    assert _select_per_replica(pooled, modes, stds, lookup) == {1: (1, 2), 2: (2, 4), 3: (3, 3)}


def test_pooled_best_is_one_of_the_per_replica_frames():
    pooled, modes, stds, lookup = _data()
    best = lookup[_select_frame(pooled, modes, stds)]
    assert best in _select_per_replica(pooled, modes, stds, lookup).values()


def test_water_penalty_is_applied_per_replica():
    pooled, modes, stds, lookup = _data()
    pen = np.zeros(len(lookup))
    pen[7] = 3.0                                   # replica 2, frame 4 lacks the water
    assert _select_per_replica(pooled, modes, stds, lookup, pen)[2] == (2, 2)


def test_pdb_parser_blank_chain_and_hybrid36_waters(tmp_path):
    from senda.analysis.distances import _parse_pdb_residues
    pdb = tmp_path / "x.pdb"
    pdb.write_text(
        "ATOM      1  N   MET     1     -30.943   9.366   2.272  1.00  0.00           N\n"
        "HETATM    2 FE1  COF   310       0.000   0.000   0.000  1.00  0.00          FE\n"
        "ATOM      3  O   WAT A000       1.000   0.000   0.000  1.00  0.00           O\n")
    assert _parse_pdb_residues(pdb) == [("A", 1, "MET"), ("A", 310, "COF")]

"""senda-complex alignment: CA atoms of amino acids in the alignment chains, no fixed residue cutoff."""

import numpy as np
import pytest

from senda.michaelis_complex.align import _parse_ca, compute_superposition


def _line(rec, name, res, ch, num, xyz):
    return f"{rec:<6}{1:5d} {name:<4} {res:>3} {ch}{num:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00\n"


def _structure(shift=(0.0, 0.0, 0.0), n=400):
    rng = np.random.default_rng(1)
    xyz = rng.normal(0, 10, (n, 3))
    lines = [_line("ATOM", "CA", "ALA", "A", i + 1, xyz[i] + shift) for i in range(n)]
    lines.append(_line("HETATM", "CA", "CA", "A", 401, (99.0, 0.0, 0.0)))            # calcium ion
    lines.append(_line("ATOM", "CA", "ALA", "B", 1, (50.0, 50.0, 50.0)))              # other chain
    return lines


def test_all_residues_of_the_alignment_chain_no_cutoff():
    keys, _ = _parse_ca(_structure())
    assert len(keys) == 400 and keys[-1] == ("A", 400)      # beyond the old 301 cutoff


def test_calcium_and_other_chains_are_excluded():
    keys, _ = _parse_ca(_structure())
    assert ("A", 401) not in keys and all(ch == "A" for ch, _ in keys)
    assert len(_parse_ca(_structure(), chains=("A", "B"))[0]) == 401


def test_residue_range_reproduces_a_fixed_window():
    keys, _ = _parse_ca(_structure(), resnum_range=(1, 301))
    assert len(keys) == 301


def test_superposition_recovers_a_translation():
    pytest.importorskip("Bio")
    rot, tran, rms = compute_superposition(_structure(shift=(1.0, 2.0, 3.0)), _structure())
    assert rms < 1e-6 and np.allclose(tran, (-1.0, -2.0, -3.0), atol=1e-6)

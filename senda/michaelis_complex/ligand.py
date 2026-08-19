"""
senda.michaelis_complex.ligand
Ligand extraction and insertion utilities for PDB lines.
"""

from __future__ import annotations


def extract_ligand(
    lines:   list[str],
    resname: str,
    chains:  set[str] | None = None,
) -> list[str]:
    """
    Return all HETATM and their paired ANISOU lines for *resname*.
    If *chains* is given, only return lines from those chains.
    """
    out: list[str] = []
    for line in lines:
        rec = line[:6]
        if rec not in ("HETATM", "ANISOU"):
            continue
        if line[17:20].strip() != resname:
            continue
        if chains is not None and line[21] not in chains:
            continue
        out.append(line)
    return out


def ligand_chains(lines: list[str], resname: str) -> set[str]:
    """Return the set of chain IDs found among HETATM lines for *resname*."""
    found: set[str] = set()
    for line in lines:
        if line[:6] == "HETATM" and line[17:20].strip() == resname:
            found.add(line[21])
    return found


def relabel_chain(lines: list[str], new_chain: str) -> list[str]:
    """Return *lines* with the chain ID column (PDB column 22) set to *new_chain*."""
    return [line[:21] + new_chain + line[22:] for line in lines]


def remove_ligand(lines: list[str], resname: str) -> list[str]:
    """Remove all HETATM and ANISOU lines whose residue name matches *resname*."""
    return [
        l for l in lines
        if not (l[:6] in ("HETATM", "ANISOU") and l[17:20].strip() == resname)
    ]


def protein_chains(
    lines:  list[str],
    chains: frozenset[str] = frozenset(("A", "B")),
) -> set[str]:
    """Return the set of chain IDs found in ATOM records, restricted to *chains*."""
    found: set[str] = set()
    for line in lines:
        if line.startswith("ATOM  "):
            ch = line[21]
            if ch in chains:
                found.add(ch)
    return found


def renumber_serial(lines: list[str]) -> list[str]:
    """
    Renumber ATOM/HETATM/ANISOU serial numbers sequentially from 1.
    Each ANISOU receives the same serial as the preceding ATOM/HETATM.
    """
    out:    list[str] = []
    serial: int       = 0

    for line in lines:
        rec = line[:6]
        if rec in ("ATOM  ", "HETATM"):
            serial += 1
            line    = rec + f"{serial:5d}" + line[11:]
        elif rec == "ANISOU":
            line    = rec + f"{serial:5d}" + line[11:]
        out.append(line)

    return out

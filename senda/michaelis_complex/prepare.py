"""
senda.michaelis_complex.prepare
Orchestrate Michaelis complex preparation for one or all mutants.

For each raw crystal-structure PDB:

  1. clean()     — strip non-A/B chains, GOL, CONECT, rename HIS, rename 216→LER
  2. align()     — Cα superposition onto the WT_LER reference; apply transform
  3. write three output files:
       {mutant}_LER_dimer.pdb  — protein + LER (from the raw structure)
       {mutant}_NIR_dimer.pdb  — protein + NIR (from the WT_NIR reference)
       {mutant}_APO_dimer.pdb  — protein only

NIR canonical coordinates are taken directly from the WT_NIR reference PDB
(which is already in the common WT_LER coordinate frame).  Each chain present
in the mutant protein receives its own copy.
"""

from __future__ import annotations

import logging
from pathlib import Path

from .clean  import clean
from .align  import compute_superposition, apply_transform
from .ligand import (
    extract_ligand,
    remove_ligand,
    protein_chains,
    renumber_serial,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _complete_missing_residues(
    aligned: list[str],
    ref:     list[str],
) -> tuple[list[str], dict[str, set[int]]]:
    """
    Copy ATOM/ANISOU records for protein residues present in *ref* but absent
    in *aligned*.

    Since *aligned* has already been superimposed onto *ref*, the reference
    coordinates are valid substitutes for structurally missing residues.

    Only protein residues (those backed by an ATOM record in *ref*) are
    considered — ligand ANISOU at e.g. residue 401 are not treated as missing
    protein residues.

    Returns the extended list and a dict {chain: {resnum, …}} of what was
    added, for logging purposes.
    """
    # Residues present in the aligned structure (protein ATOM only)
    present: set[tuple[str, int]] = set()
    for line in aligned:
        if line.startswith("ATOM  ") and line[21] in ("A", "B"):
            try:
                present.add((line[21], int(line[22:26])))
            except ValueError:
                pass

    # Protein residues available in the reference (ATOM records only)
    ref_protein: set[tuple[str, int]] = set()
    for line in ref:
        if line.startswith("ATOM  ") and line[21] in ("A", "B"):
            try:
                ref_protein.add((line[21], int(line[22:26])))
            except ValueError:
                pass

    # Residues to fill = in ref protein but missing from aligned
    to_fill = ref_protein - present

    # Collect ATOM and ANISOU lines for those residues
    added:      list[str]           = []
    added_keys: dict[str, set[int]] = {"A": set(), "B": set()}

    for line in ref:
        if line[:6] not in ("ATOM  ", "ANISOU"):
            continue
        ch = line[21]
        if ch not in ("A", "B"):
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if (ch, resnum) in to_fill:
            added.append(line)
            added_keys[ch].add(resnum)

    return aligned + added, added_keys


def _sort_protein_lines(lines: list[str]) -> list[str]:
    """
    Sort ATOM/ANISOU lines by residue number.

    Atoms within the same residue keep their original relative order
    (N, CA, C, O, CB, …) via a stable sort on (resnum, original_index).
    This is needed when missing residues copied from the reference are
    appended at the end and must be interleaved with the existing ones.
    """
    keyed = []
    for i, line in enumerate(lines):
        try:
            resnum = int(line[22:26])
        except ValueError:
            resnum = 0
        keyed.append((resnum, i, line))
    keyed.sort(key=lambda x: (x[0], x[1]))
    return [line for _, _, line in keyed]


def _write_pdb(
    body:         list[str],
    ligand_lines: list[str],
    out_path:     Path,
) -> None:
    """
    Write a PDB file with TER records between protein chains.

    Output order:
      <header records>
      <chain A ATOM/ANISOU — sorted by residue number>
      TER
      <chain B ATOM/ANISOU — sorted by residue number>   (omitted when absent)
      TER
      <ligand HETATM/ANISOU>
      <water HETATM/ANISOU>
      TER

    Sorting by residue number is required when missing residues have been
    copied from the reference and appended at the end of the body.

    tleap requires TER between chains; without it it tries to bond the
    C-terminus of chain A to the N-terminus of chain B.
    """
    header: list[str]            = []
    prot:   dict[str, list[str]] = {"A": [], "B": []}
    water:  list[str]            = []

    for line in body:
        if line.startswith(("TER", "END")):
            continue
        rec = line[:6]
        if rec in ("ATOM  ", "ANISOU"):
            ch = line[21]
            if ch in prot:
                prot[ch].append(line)
        elif rec == "HETATM":
            if line[17:20].strip() == "HOH":
                water.append(line)
        else:
            header.append(line)

    out: list[str] = list(header)
    for ch in ("A", "B"):
        if prot[ch]:
            out.extend(_sort_protein_lines(prot[ch]))
            out.append("TER\n")
    out.extend(ligand_lines)
    out.extend(water)
    out.append("TER\n")

    out_path.write_text("".join(renumber_serial(out)))


# ---------------------------------------------------------------------------
# Single-structure preparation
# ---------------------------------------------------------------------------

def prepare(
    raw_pdb:    Path,
    ref_ler:    list[str],
    ref_nir:    list[str],
    output_dir: Path,
    mutant:     str,
    his_rename: dict[int, str] | None = None,
) -> None:
    """
    Produce *{mutant}_LER_dimer.pdb*, *{mutant}_NIR_dimer.pdb*, and
    *{mutant}_APO_dimer.pdb* in *output_dir* from a raw crystal-structure PDB.

    Parameters
    ----------
    raw_pdb    : raw crystal structure (e.g. wildtype_1216_refmac6.pdb)
    ref_ler    : lines of WT_LER_dimer.pdb (pre-loaded) used as superposition reference
    ref_nir    : lines of WT_NIR_dimer.pdb (pre-loaded) — source of NIR coordinates
    output_dir : directory where the three output PDBs are written
    mutant     : name prefix for output files (e.g. "WT", "E166V")
    his_rename : residue-number → new name mapping (defaults: Mpro protonation)
    """
    raw_lines = raw_pdb.read_text().splitlines(keepends=True)

    # 1. Clean
    cleaned = clean(raw_lines, his_rename)

    # 2. Align to WT_LER reference frame
    rot, tran, rms = compute_superposition(cleaned, ref_ler)
    log.info("[%s] Cα RMSD to reference: %.3f Å", mutant, rms)
    aligned = apply_transform(cleaned, rot, tran)

    # 3. Fill residues missing from the crystal structure using the reference
    aligned, filled = _complete_missing_residues(aligned, ref_ler)
    for ch, resnums in filled.items():
        if resnums:
            log.info("[%s] chain %s: copied %d missing residues from reference (%s)",
                     mutant, ch, len(resnums), sorted(resnums))
            print(f"  chain {ch}: filled residues {sorted(resnums)} from reference")

    chains = protein_chains(aligned)

    # 3a. LER complex — LER is already in aligned lines (renamed from 216)
    ler_ligand   = extract_ligand(aligned, "LER")
    protein_body = remove_ligand(aligned, "LER")
    output_dir.mkdir(parents=True, exist_ok=True)
    _write_pdb(protein_body, ler_ligand, output_dir / f"{mutant}_LER_dimer.pdb")

    # 3b. NIR complex — canonical NIR from WT_NIR reference per available chain
    nir_ligand: list[str] = []
    for ch in sorted(chains):
        nir_ligand.extend(extract_ligand(ref_nir, "NIR", chains={ch}))

    if not nir_ligand:
        log.warning("[%s] No NIR found in reference for chains %s", mutant, chains)

    _write_pdb(protein_body, nir_ligand, output_dir / f"{mutant}_NIR_dimer.pdb")

    # 3c. APO — protein only
    _write_pdb(protein_body, [], output_dir / f"{mutant}_APO_dimer.pdb")


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def prepare_all(
    structures:   dict[str, str],
    raw_pdbs_dir: Path,
    ref_ler_pdb:  Path,
    ref_nir_pdb:  Path,
    output_dir:   Path,
    his_rename:   dict[int, str] | None = None,
    force:        bool = False,
) -> None:
    """
    Run :func:`prepare` for each entry in *structures*.

    Parameters
    ----------
    structures   : mapping of {raw_pdb_filename: mutant_name}
    raw_pdbs_dir : directory containing the raw PDB files
    ref_ler_pdb  : WT_LER_dimer.pdb (alignment reference; read once before any writes)
    ref_nir_pdb  : WT_NIR_dimer.pdb (NIR coordinate source)
    output_dir   : where to write the curated dimers
    his_rename   : HIS protonation overrides
    force        : overwrite existing output files
    """
    # Read reference files once upfront so that overwriting WT outputs
    # (when WT is also listed in structures) doesn't corrupt the reference data.
    ref_ler = ref_ler_pdb.read_text().splitlines(keepends=True)
    ref_nir = ref_nir_pdb.read_text().splitlines(keepends=True)

    for filename, mutant in structures.items():
        raw_pdb = raw_pdbs_dir / filename
        if not raw_pdb.exists():
            raise FileNotFoundError(f"Raw PDB not found: {raw_pdb}")

        outputs = [
            output_dir / f"{mutant}_LER_dimer.pdb",
            output_dir / f"{mutant}_NIR_dimer.pdb",
            output_dir / f"{mutant}_APO_dimer.pdb",
        ]
        if not force and all(p.exists() for p in outputs):
            log.info("[%s] already done — skipping (use --force to redo)", mutant)
            print(f"[{mutant}] already done — skipping (use --force to redo)")
            continue

        print(f"[{mutant}] {raw_pdb.name} ...")
        prepare(raw_pdb, ref_ler, ref_nir, output_dir, mutant, his_rename)
        print(f"[{mutant}] done → {output_dir}/")

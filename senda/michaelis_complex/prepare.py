"""
senda.michaelis_complex.prepare
Orchestrate Michaelis complex preparation for one or all mutants.

For each raw crystal-structure PDB:

  1. clean()     — strip unwanted chains/residues, apply residue_renames,
                   apply protonation states from reference structure
  2. align()     — Cα superposition onto the alignment reference; apply transform
  3. write output files:
       {mutant}_{inhibitor}_dimer.pdb  — one per inhibitor in inhibitor_sources
       {mutant}_APO_dimer.pdb          — protein only (always written)

Inhibitor coordinates come from two sources, configured per inhibitor:
  "native"          — the ligand is already present in the raw PDB (after
                      residue_renames are applied)
  "/path/to/ref.pdb" — coordinates are copied from the named reference PDB,
                       one copy per chain present in the mutant protein
"""

from __future__ import annotations

import logging
from pathlib import Path

from .clean  import clean, build_protonation_map
from .align  import compute_superposition, apply_transform
from .ligand import (
    extract_ligand,
    remove_ligand,
    protein_chains,
    ligand_chains,
    relabel_chain,
    renumber_serial,
)

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _complete_missing_residues(
    aligned: list[str],
    ref:     list[str],
    chains:  frozenset[str],
) -> tuple[list[str], dict[str, set[int]]]:
    """
    Copy ATOM/ANISOU records for protein residues present in *ref* but absent
    in *aligned*.

    Since *aligned* has already been superimposed onto *ref*, the reference
    coordinates are valid substitutes for structurally missing residues.
    """
    present: set[tuple[str, int]] = set()
    for line in aligned:
        if line.startswith("ATOM  ") and line[21] in chains:
            try:
                present.add((line[21], int(line[22:26])))
            except ValueError:
                pass

    ref_protein: set[tuple[str, int]] = set()
    for line in ref:
        if line.startswith("ATOM  ") and line[21] in chains:
            try:
                ref_protein.add((line[21], int(line[22:26])))
            except ValueError:
                pass

    to_fill = ref_protein - present

    added:      list[str]           = []
    added_keys: dict[str, set[int]] = {ch: set() for ch in chains}

    for line in ref:
        if line[:6] not in ("ATOM  ", "ANISOU"):
            continue
        ch = line[21]
        if ch not in chains:
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        if (ch, resnum) in to_fill:
            added.append(line)
            if ch in added_keys:
                added_keys[ch].add(resnum)

    return aligned + added, added_keys


def _trim_to_reference(
    lines:  list[str],
    ref:    list[str],
    chains: frozenset[str],
    mutant: str,
) -> list[str]:
    """
    Remove ATOM records whose residue number falls outside the per-chain
    range defined by the alignment reference.  HETATM records (ligands,
    crystal waters) are left untouched.
    """
    ref_range: dict[str, tuple[int, int]] = {}
    for line in ref:
        if not line.startswith("ATOM  "):
            continue
        ch = line[21]
        if ch not in chains:
            continue
        try:
            resnum = int(line[22:26])
        except ValueError:
            continue
        lo, hi = ref_range.get(ch, (resnum, resnum))
        ref_range[ch] = (min(lo, resnum), max(hi, resnum))

    result   = []
    trimmed: dict[str, set[int]] = {ch: set() for ch in chains}
    for line in lines:
        if line.startswith("ATOM  "):
            ch = line[21]
            if ch in ref_range:
                try:
                    resnum = int(line[22:26])
                except ValueError:
                    result.append(line)
                    continue
                lo, hi = ref_range[ch]
                if resnum < lo or resnum > hi:
                    trimmed[ch].add(resnum)
                    continue
        result.append(line)

    for ch, resnums in trimmed.items():
        if resnums:
            log.info("[%s] chain %s: trimmed %d residues outside reference range (%s)",
                     mutant, ch, len(resnums), sorted(resnums))
            print(f"  chain {ch}: trimmed residues {sorted(resnums)} (outside reference range)")

    return result


def _sort_protein_lines(lines: list[str]) -> list[str]:
    """
    Sort ATOM/ANISOU lines by residue number, preserving within-residue order.
    Needed after missing residues copied from the reference are appended at the end.
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
    chains:       frozenset[str],
) -> None:
    """
    Write a PDB file ordered as:
      protein chain A  ->  TER
      protein chain B  ->  TER
      inhibitor chain A  ->  TER
      inhibitor chain B  ->  TER
      water  ->  TER

    ANISOU records are dropped throughout (protein, ligand, and water).
    tleap requires TER between chains and between ligand copies.
    """
    header: list[str]            = []
    prot:   dict[str, list[str]] = {ch: [] for ch in chains}
    water:  list[str]            = []

    for line in body:
        if line.startswith(("TER", "END", "ANISOU")):
            continue
        rec = line[:6]
        if rec == "ATOM  ":
            ch = line[21]
            if ch in prot:
                prot[ch].append(line)
        elif rec == "HETATM":
            if line[17:20].strip() == "HOH":
                water.append(line)
        else:
            header.append(line)

    out: list[str] = list(header)
    for ch in sorted(chains):
        if prot[ch]:
            out.extend(_sort_protein_lines(prot[ch]))
            out.append("TER\n")

    prev_ch = None
    for line in ligand_lines:
        if line[:6] == "ANISOU":
            continue
        if line[:6] == "HETATM":
            ch = line[21]
            if prev_ch is not None and ch != prev_ch:
                out.append("TER\n")
            prev_ch = ch
        out.append(line)
    if prev_ch is not None:
        out.append("TER\n")

    out.extend(water)
    out.append("TER\n")

    out_path.write_text("".join(renumber_serial(out)))


# ---------------------------------------------------------------------------
# Single-structure preparation
# ---------------------------------------------------------------------------

def prepare(
    raw_pdb:           Path,
    alignment_ref:     list[str],
    inhibitor_sources: dict[str, str],
    ref_inh_lines:     dict[str, list[str]],
    output_dir:        Path,
    mutant:            str,
    prot_map:          dict[tuple[str, int], str] | None = None,
    residue_renames:   list[tuple[str, str]] | None = None,
    chains:            frozenset[str] = frozenset(("A", "B")),
    include_apo:       bool = True,
    disulfides:        list[tuple[str, int, str, int]] | None = None,
) -> None:
    """
    Produce inhibitor complex and APO PDBs in *output_dir*.

    Parameters
    ----------
    raw_pdb           : raw crystal structure to process
    alignment_ref     : lines of the alignment reference PDB (pre-loaded)
    inhibitor_sources : {inhibitor_name: "native" | path_string}
    ref_inh_lines     : pre-loaded lines for each non-native inhibitor reference
    output_dir        : directory where output PDBs are written
    mutant            : name prefix for output files (e.g. "WT", "E166V")
    prot_map          : (chain, resnum) → resname from reference structure
    residue_renames   : [(from, to), …] applied before protonation
    chains            : chain IDs to retain
    disulfides        : [(chain1, resnum1, chain2, resnum2), …] confirmed
                        disulfide pairs, renamed to CYX
    """
    raw_lines = raw_pdb.read_text().splitlines(keepends=True)

    # 1. Clean: renames + protonation
    cleaned = clean(raw_lines, prot_map, residue_renames, chains, disulfides=disulfides)

    # 2. Align to reference frame
    rot, tran, rms = compute_superposition(cleaned, alignment_ref)
    log.info("[%s] Cα RMSD to reference: %.3f Å", mutant, rms)
    aligned = apply_transform(cleaned, rot, tran)

    # 3. Fill residues missing from the crystal structure using the reference
    aligned, filled = _complete_missing_residues(aligned, alignment_ref, chains)
    for ch, resnums in filled.items():
        if resnums:
            log.info("[%s] chain %s: copied %d missing residues from reference (%s)",
                     mutant, ch, len(resnums), sorted(resnums))
            print(f"  chain {ch}: filled residues {sorted(resnums)} from reference")

    # 3b. Trim residues outside the reference range (some crystal structures
    #     resolve more C-terminal residues than the reference, which would
    #     introduce spurious charge differences between systems).
    aligned = _trim_to_reference(aligned, alignment_ref, chains, mutant)

    chains_present = protein_chains(aligned, chains)

    # Protein body = aligned lines with all native ligands removed
    native_names = [n for n, s in inhibitor_sources.items() if s == "native"]
    protein_body = aligned
    for name in native_names:
        protein_body = remove_ligand(protein_body, name)

    output_dir.mkdir(parents=True, exist_ok=True)

    # 4. Write one PDB per inhibitor
    for inh_name, source in inhibitor_sources.items():
        if source == "native":
            ligand_lines = extract_ligand(aligned, inh_name)
        else:
            # Match target chains to the reference's ligand chains
            # positionally (1st retained chain <-> 1st ligand instance in
            # the reference, 2nd <-> 2nd, ...) and relabel the copied lines
            # to the target chain letter. The reference PDB isn't required
            # to use the same chain letters as the target -- e.g. a
            # reference with ligand chains {A, B} still supplies target
            # chains {A, C} correctly.
            ref_lines     = ref_inh_lines[inh_name]
            target_chains = sorted(chains_present)
            ref_chains    = sorted(ligand_chains(ref_lines, inh_name))
            if len(ref_chains) < len(target_chains):
                log.warning(
                    "[%s] Reference for %s has ligand in %d chain(s) %s but "
                    "%d target chain(s) %s need one -- some target chains "
                    "will have no %s copy.",
                    mutant, inh_name, len(ref_chains), ref_chains,
                    len(target_chains), target_chains, inh_name,
                )
            ligand_lines = []
            for target_ch, ref_ch in zip(target_chains, ref_chains):
                matched = extract_ligand(ref_lines, inh_name, chains={ref_ch})
                ligand_lines.extend(relabel_chain(matched, target_ch))
            if not ligand_lines:
                log.warning("[%s] No %s found in reference for chains %s",
                            mutant, inh_name, chains_present)

        _write_pdb(protein_body, ligand_lines,
                   output_dir / f"{mutant}_{inh_name}_dimer.pdb", chains)

    # 5. APO — protein only (only when requested)
    if include_apo:
        _write_pdb(protein_body, [], output_dir / f"{mutant}_APO_dimer.pdb", chains)


# ---------------------------------------------------------------------------
# Batch runner
# ---------------------------------------------------------------------------

def prepare_all(
    structures:           dict[str, str],
    raw_pdbs_dir:         Path,
    alignment_ref_pdb:    Path,
    inhibitor_sources:    dict[str, str],
    reference_enzyme_pdb: Path | None,
    output_dir:           Path,
    chains:               frozenset[str] = frozenset(("A", "B")),
    residue_renames:      list[tuple[str, str]] | None = None,
    force:                bool = False,
    include_apo:          bool = True,
    disulfides:           list[tuple[str, int, str, int]] | None = None,
) -> None:
    """
    Run :func:`prepare` for each entry in *structures*.

    Parameters
    ----------
    structures            : {raw_pdb_filename: mutant_name}
    raw_pdbs_dir          : directory containing the raw PDB files
    alignment_ref_pdb     : PDB used for Cα superposition and missing-residue filling
    inhibitor_sources     : {inhibitor_name: "native" | path_to_ref_pdb}
    reference_enzyme_pdb  : PDB from which protonation states are read; ligands
                            and water are ignored automatically so a holo structure
                            is fine. None → no protonation renaming.
    output_dir            : where to write the curated PDBs
    chains                : chain IDs to retain
    residue_renames       : [(from, to), …] applied before protonation
    force                 : overwrite existing output files
    disulfides            : [(chain1, resnum1, chain2, resnum2), …] confirmed
                            disulfide pairs, renamed to CYX
    """
    # Read alignment reference once upfront so that overwriting WT outputs
    # does not corrupt it when WT is also listed in structures.
    alignment_ref = alignment_ref_pdb.read_text().splitlines(keepends=True)

    # Build protonation map from reference enzyme PDB (protein ATOM records only)
    prot_map: dict[tuple[str, int], str] | None = None
    if reference_enzyme_pdb is not None:
        ref_enz = reference_enzyme_pdb.read_text().splitlines(keepends=True)
        prot_map = build_protonation_map(ref_enz, chains)
        log.info("Protonation map loaded from %s (%d positions)",
                 reference_enzyme_pdb.name, len(prot_map))

    # Pre-load reference PDBs for non-native inhibitors
    ref_inh_lines: dict[str, list[str]] = {}
    for inh_name, source in inhibitor_sources.items():
        if source != "native":
            ref_inh_lines[inh_name] = Path(source).read_text().splitlines(keepends=True)

    # Expected outputs per structure
    inh_names = list(inhibitor_sources.keys())

    for filename, mutant in structures.items():
        raw_pdb = raw_pdbs_dir / filename
        if not raw_pdb.exists():
            raise FileNotFoundError(f"Raw PDB not found: {raw_pdb}")

        outputs = [output_dir / f"{mutant}_{n}_dimer.pdb" for n in inh_names]
        if include_apo:
            outputs.append(output_dir / f"{mutant}_APO_dimer.pdb")
        if not force and all(p.exists() for p in outputs):
            log.info("[%s] already done — skipping (use --force to redo)", mutant)
            print(f"[{mutant}] already done — skipping (use --force to redo)")
            continue

        print(f"[{mutant}] {raw_pdb.name} ...")
        prepare(
            raw_pdb, alignment_ref,
            inhibitor_sources, ref_inh_lines,
            output_dir, mutant,
            prot_map, residue_renames, chains,
            include_apo=include_apo,
            disulfides=disulfides,
        )
        print(f"[{mutant}] done → {output_dir}/")

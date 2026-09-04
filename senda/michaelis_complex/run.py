"""
senda.michaelis_complex.run
CLI entry point for senda-complex.

Usage:
    senda-complex --config config.yaml
    senda-complex --config config.yaml --force
    senda-complex --config config.yaml --report-disulfides
"""

import argparse
import logging
import sys
from collections import Counter
from pathlib import Path

log = logging.getLogger(__name__)


def _validate_disulfides(disulfides: list[tuple[str, int, str, int]]) -> None:
    """Hard-stop if any (chain, resnum) endpoint appears in more than one pair."""
    counts = Counter()
    for c1, r1, c2, r2 in disulfides:
        counts[(c1, r1)] += 1
        counts[(c2, r2)] += 1
    repeated = [key for key, n in counts.items() if n > 1]
    if repeated:
        sys.exit(
            "Error: michaelis_complex.disulfides lists the same residue in more "
            f"than one pair: {repeated} — a cysteine can only form one disulfide bond."
        )


def main() -> None:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
    )

    parser = argparse.ArgumentParser(
        description=(
            "Prepare Michaelis complexes from raw crystal structures.\n"
            "Reads michaelis_complex: section from config.yaml."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", required=True, metavar="CONFIG",
                        help="Path to senda config.yaml")
    parser.add_argument("--force", action="store_true",
                        help="Re-run even if output files already exist")
    parser.add_argument("--report-disulfides", action="store_true",
                        help="Detect candidate Cys SG-SG disulfide pairs per raw "
                             "structure and exit (no output files written)")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        sys.exit(f"Error: config file not found: {config_path}")

    try:
        import yaml
    except ImportError:
        sys.exit("pyyaml is required: pip install pyyaml")

    with open(config_path) as fh:
        raw = yaml.safe_load(fh) or {}

    mc = raw.get("michaelis_complex")
    if not mc:
        sys.exit("Error: missing 'michaelis_complex:' section in config")

    cwd = Path.cwd()

    raw_pdbs_dir = cwd / mc.get("raw_pdbs_dir", "raw_pdbs")
    output_dir   = cwd / mc.get("output_dir",   "protein")
    structures = mc.get("structures", {})

    if not structures:
        sys.exit("Error: no structures listed under michaelis_complex.structures")

    # Filter to only the mutants listed at the top level, if specified.
    mutants = raw.get("mutants") or []
    if mutants:
        structures = {f: m for f, m in structures.items() if m in mutants}
        if not structures:
            sys.exit(
                "Error: none of the mutants listed under 'mutants' match any entry "
                "in michaelis_complex.structures"
            )

    # alignment_reference is required
    if "alignment_reference" not in mc:
        sys.exit("Error: michaelis_complex.alignment_reference is required")
    alignment_ref_pdb = cwd / mc["alignment_reference"]

    # reference_enzyme_pdb defaults to alignment_reference
    ref_enzyme_key    = mc.get("reference_enzyme_pdb", mc["alignment_reference"])
    reference_enzyme_pdb = cwd / ref_enzyme_key

    chains = frozenset(str(c) for c in mc.get("chains", ["A", "B"]))

    if args.report_disulfides:
        from .disulfides import find_disulfide_candidates

        all_candidates: list[dict] = []
        for filename, mutant in structures.items():
            raw_pdb = raw_pdbs_dir / filename
            if not raw_pdb.exists():
                sys.exit(f"Error: not found: {raw_pdb}")
            lines = raw_pdb.read_text().splitlines(keepends=True)
            candidates = find_disulfide_candidates(lines, chains)

            print(f"[{mutant}] {filename}")
            if not candidates:
                print("  (no candidates within 3.0 A)")
            else:
                print(f"  {'Chain1':<7}{'Res1':<6}{'Chain2':<7}{'Res2':<6}{'Distance':<10}Band")
                for c in candidates:
                    flag = "  <- verify manually" if c["band"] == "borderline" else ""
                    print(f"  {c['chain1']:<7}{c['resnum1']:<6}{c['chain2']:<7}"
                          f"{c['resnum2']:<6}{c['distance']:<9.2f} A {c['band']}{flag}")
                all_candidates.extend(c for c in candidates if c["band"] == "bonded")
            print()

        if all_candidates:
            seen = set()
            print("Add confirmed pairs to michaelis_complex.disulfides in config.yaml:")
            print("  michaelis_complex:")
            print("    disulfides:")
            for c in all_candidates:
                key = (c["chain1"], c["resnum1"], c["chain2"], c["resnum2"])
                if key in seen:
                    continue
                seen.add(key)
                print(f"      - chain1: {c['chain1']}")
                print(f"        resnum1: {c['resnum1']}")
                print(f"        chain2: {c['chain2']}")
                print(f"        resnum2: {c['resnum2']}")
        sys.exit(0)

    residue_renames = [
        (str(r["from"]), str(r["to"]))
        for r in mc.get("residue_renames", [])
    ]

    disulfides = [
        (str(d["chain1"]), int(d["resnum1"]), str(d["chain2"]), int(d["resnum2"]))
        for d in mc.get("disulfides", [])
    ]
    _validate_disulfides(disulfides)

    # inhibitor_sources: {name: "native"} or {name: path_relative_to_cwd}
    raw_sources = mc.get("inhibitor_sources", {})
    if not raw_sources:
        sys.exit("Error: michaelis_complex.inhibitor_sources is required")
    inhibitor_sources: dict[str, str] = {}
    for name, source in raw_sources.items():
        if str(source).lower() == "native":
            inhibitor_sources[name] = "native"
        else:
            inhibitor_sources[name] = str(cwd / source)

    # Filter to only the inhibitors listed at the top level, if specified.
    # APO is handled separately (include_apo flag) and has no inhibitor_sources entry.
    inhibitors  = raw.get("inhibitors") or []
    include_apo = not inhibitors or "APO" in inhibitors
    if inhibitors:
        ligand_inh        = [i for i in inhibitors if i != "APO"]
        inhibitor_sources = {n: s for n, s in inhibitor_sources.items() if n in ligand_inh}
        if ligand_inh and not inhibitor_sources:
            sys.exit(
                "Error: none of the inhibitors listed under 'inhibitors' match any entry "
                "in michaelis_complex.inhibitor_sources"
            )

    # Validate paths
    for path in (raw_pdbs_dir, alignment_ref_pdb, reference_enzyme_pdb):
        if not path.exists():
            sys.exit(f"Error: not found: {path}")
    for name, source in inhibitor_sources.items():
        if source != "native" and not Path(source).exists():
            sys.exit(f"Error: inhibitor reference not found for {name}: {source}")

    print(f"Raw PDBs dir      : {raw_pdbs_dir}")
    print(f"Output dir        : {output_dir}")
    print(f"Alignment ref     : {alignment_ref_pdb}")
    print(f"Reference enzyme  : {reference_enzyme_pdb}")
    print(f"Chains            : {sorted(chains)}")
    print(f"Residue renames   : {residue_renames or '(none)'}")
    print(f"Disulfides        : {disulfides or '(none)'}")
    print(f"Inhibitor sources : {inhibitor_sources}")
    print(f"Structures        : {len(structures)}")
    print()

    from .prepare import prepare_all

    try:
        prepare_all(
            structures=structures,
            raw_pdbs_dir=raw_pdbs_dir,
            alignment_ref_pdb=alignment_ref_pdb,
            inhibitor_sources=inhibitor_sources,
            reference_enzyme_pdb=reference_enzyme_pdb,
            output_dir=output_dir,
            chains=chains,
            residue_renames=residue_renames,
            force=args.force,
            include_apo=include_apo,
            disulfides=disulfides,
        )
    except FileNotFoundError as exc:
        sys.exit(f"Error: {exc}")
    except Exception as exc:
        log.exception("Unexpected error")
        sys.exit(f"Error: {exc}")

    print(f"\nAll done. Outputs in: {output_dir}")

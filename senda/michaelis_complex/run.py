"""
senda.michaelis_complex.run
CLI entry point for senda-complex.

Usage:
    senda-complex --config config.yaml
    senda-complex --config config.yaml --force
"""

import argparse
import logging
import sys
from pathlib import Path

log = logging.getLogger(__name__)


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

    residue_renames = [
        (str(r["from"]), str(r["to"]))
        for r in mc.get("residue_renames", [])
    ]

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
        )
    except FileNotFoundError as exc:
        sys.exit(f"Error: {exc}")
    except Exception as exc:
        log.exception("Unexpected error")
        sys.exit(f"Error: {exc}")

    print(f"\nAll done. Outputs in: {output_dir}")

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

_DEFAULT_HIS: dict[int, str] = {
    41:  "HID",
    64:  "HIE",
    80:  "HID",
    163: "HIE",
    164: "HIE",
    172: "HIE",
    246: "HIE",
}


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

    raw_pdbs_dir  = cwd / mc.get("raw_pdbs_dir",  "raw_pdbs")
    output_dir    = cwd / mc.get("output_dir",     "protein")
    structures    = mc.get("structures", {})

    try:
        ref_ler_pdb  = cwd / mc["reference_ler"]
        ref_nir_pdb  = cwd / mc["reference_nir"]
    except KeyError as exc:
        sys.exit(f"Error: missing required key in michaelis_complex config: {exc}")

    if not structures:
        sys.exit("Error: no structures listed under michaelis_complex.structures")

    # Build HIS rename map: start from defaults, apply config overrides
    his_cfg    = mc.get("his_rename", {}) or {}
    his_rename = dict(_DEFAULT_HIS)
    his_rename.update({int(k): str(v) for k, v in his_cfg.items()})

    for path in (raw_pdbs_dir, ref_ler_pdb, ref_nir_pdb):
        if not path.exists():
            sys.exit(f"Error: not found: {path}")

    print(f"Raw PDBs dir  : {raw_pdbs_dir}")
    print(f"Output dir    : {output_dir}")
    print(f"Reference LER : {ref_ler_pdb}")
    print(f"Reference NIR : {ref_nir_pdb}")
    print(f"Structures    : {len(structures)}")
    print()

    from .prepare import prepare_all

    try:
        prepare_all(
            structures=structures,
            raw_pdbs_dir=raw_pdbs_dir,
            ref_ler_pdb=ref_ler_pdb,
            ref_nir_pdb=ref_nir_pdb,
            output_dir=output_dir,
            his_rename=his_rename,
            force=args.force,
        )
    except FileNotFoundError as exc:
        sys.exit(f"Error: {exc}")
    except Exception as exc:
        log.exception("Unexpected error")
        sys.exit(f"Error: {exc}")

    print(f"\nAll done. Outputs in: {output_dir}")

import sys
from datetime import datetime
from pathlib import Path


class _Tee:
    """Write to two streams simultaneously (terminal + log file)."""
    def __init__(self, primary, secondary):
        self._primary   = primary
        self._secondary = secondary

    def write(self, data):
        self._primary.write(data)
        self._secondary.write(data)

    def flush(self):
        self._primary.flush()
        self._secondary.flush()

    # Delegate attribute lookups (e.g. .encoding) to the primary stream
    def __getattr__(self, name):
        return getattr(self._primary, name)


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Analyse reactive distances from NVT trajectories and select "
                    "a representative frame per chain for QM/MM input."
    )
    parser.add_argument("--config", required=True, metavar="CONFIG",
                        help="Path to senda config.yaml")
    parser.add_argument("--inhibitors", nargs="+", metavar="INH",
                        help="Inhibitors to analyse (default: all from config)")
    parser.add_argument("--mutants", nargs="+", metavar="MUT",
                        help="Mutants to analyse (default: all from config)")
    args = parser.parse_args()

    config_path = Path(args.config)
    if not config_path.exists():
        sys.exit(f"Error: config not found: {config_path}")

    try:
        import yaml
    except ImportError:
        sys.exit("pyyaml is required: pip install pyyaml")

    with open(config_path) as fh:
        cfg = yaml.safe_load(fh) or {}

    analysis_cfg = cfg.get("analysis") or {}

    inhibitors = args.inhibitors or cfg.get("inhibitors") or []
    mutants    = args.mutants    or cfg.get("mutants")    or []
    n_replicas = int(cfg.get("replicas", 1))
    cwd        = Path.cwd()

    mc_cfg      = cfg.get("michaelis_complex") or {}
    chains      = mc_cfg.get("chains") or ["A"]
    protein_dir = cwd / mc_cfg.get("output_dir", "protein")

    if not inhibitors:
        sys.exit("Error: no inhibitors specified (config or --inhibitors)")
    if not mutants:
        sys.exit("Error: no mutants specified (config or --mutants)")

    # Set up log file -- tee stdout so terminal output is preserved
    log_dir  = cwd / "Analysis"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / "senda_analyse.log"
    log_file = open(log_path, "w")
    log_file.write(f"senda-analyse  {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
    log_file.write(f"config: {config_path.resolve()}\n")
    log_file.write("-" * 60 + "\n")
    log_file.flush()

    original_stdout = sys.stdout
    sys.stdout = _Tee(original_stdout, log_file)

    try:
        from .distances import analyse

        for inh in inhibitors:
            if inh.upper() == "APO":
                print(f"Skipping APO (reactive distances require a ligand)")
                continue

            # Per-ligand distances take precedence; fall back to top-level list (legacy)
            inh_cfg    = analysis_cfg.get(inh) or {}
            dist_specs = inh_cfg.get("reactive_distances") \
                         or analysis_cfg.get("reactive_distances") \
                         or []
            rdf_specs              = inh_cfg.get("water_rdf") or []
            exclude_water          = inh_cfg.get("exclude_water") or []
            exclude_water_neighbor = inh_cfg.get("exclude_water_neighbor")
            if not dist_specs:
                print(f"Skipping {inh}: no reactive_distances configured")
                continue

            for mut in mutants:
                print(f"\n=== {inh} / {mut} ===")
                try:
                    analyse(inh, mut, dist_specs, rdf_specs, n_replicas, chains, protein_dir, cwd,
                            exclude_water=exclude_water,
                            exclude_water_neighbor=exclude_water_neighbor)
                except FileNotFoundError as exc:
                    print(f"  WARNING: {exc} -- skipping")
                except Exception as exc:
                    print(f"  ERROR: {exc}")
                    raise
    finally:
        sys.stdout = original_stdout
        log_file.close()
        print(f"Log saved: {log_path}")

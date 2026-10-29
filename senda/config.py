"""
senda.config

Shared config post-processing used by every senda-* entry point that reads
the slurm: section (senda-sim, senda-slurm, senda-qmmm).
"""
from __future__ import annotations


def resolve_slurm_profile(cfg: dict) -> dict:
    """
    If slurm.profiles is set, replace cfg['slurm'] with the selected
    profile's settings -- a flat dict shaped exactly like the plain
    slurm: block every senda-* command already expects, so nothing
    downstream needs to know profiles exist.

    slurm.active_profile picks which profile is in effect:

        slurm:
          active_profile: marenostrum
          profiles:
            bluepebble: {account: ..., gpu: {...}, ...}
            marenostrum: {account: ..., cpu: {...}, env_setup: ..., ...}

    Configs that don't use slurm.profiles are returned unchanged (the
    plain flat slurm: block keeps working exactly as before).
    """
    slurm = cfg.get("slurm")
    if not slurm or "profiles" not in slurm:
        return cfg

    profiles = slurm["profiles"]
    active   = slurm.get("active_profile")
    if not active:
        raise ValueError(
            "slurm.profiles is set but slurm.active_profile is not -- "
            f"choose one of: {list(profiles)}"
        )
    if active not in profiles:
        raise ValueError(
            f"slurm.active_profile {active!r} not found in slurm.profiles "
            f"(available: {list(profiles)})"
        )

    cfg = dict(cfg)
    cfg["slurm"] = profiles[active]
    return cfg


# Fields placed explicitly by templates (as __NTASKS__, __TIME__, etc. --
# QM/MM's string stage computes NTASKS_STRING from ntasks, for instance)
# rather than emitted generically -- skipped by default so they're never
# written out twice. Callers with no such reserved fields (e.g. classical
# MD, where every field is a plain passthrough) can override this with
# reserved=set().
DEFAULT_SBATCH_RESERVED = {"ntasks", "time", "time_string", "cpus-per-task"}

# Never a real SLURM flag -- purely a senda-internal key (the QM/MM
# string stage's wall-time override). Always excluded regardless of what
# a caller passes as reserved, since sbatch would reject an unrecognized
# --time_string= option outright.
_NEVER_A_FLAG = {"time_string"}


def sbatch_lines(section: dict, reserved: "set | None" = None, **extra) -> str:
    """
    Build '#SBATCH --key=value' lines for whichever fields are present in
    *section* (plus any keyword overrides), skipping *reserved* keys that
    the caller places explicitly. Missing/empty fields are simply omitted
    -- no placeholder comments, and no field is ever required.

    This is the one place cluster-specific SLURM fields get turned into
    #SBATCH lines -- don't hardcode a fixed required/optional field list
    anywhere else; different clusters need different subsets (account,
    partition, qos, gres, mem, ... are all genuinely optional depending
    on the cluster).
    """
    if reserved is None:
        reserved = DEFAULT_SBATCH_RESERVED
    excluded = reserved | _NEVER_A_FLAG
    fields = {**section, **extra}
    lines = [
        f"#SBATCH --{key}={value}"
        for key, value in fields.items()
        if key not in excluded and value not in (None, "")
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Complex PDB naming and chain selection (monomers and multimers alike)
# ---------------------------------------------------------------------------

# Amino-acid residue names (with AMBER protonation variants): a chain holding
# any of these is a protein chain. Ligands and cofactors may carry their own
# chain IDs and are not counted.
AMINO_ACIDS = frozenset({
    "ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "CYM", "GLN", "GLU",
    "GLH", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN",
    "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
})


def complex_pdb_name(mut: str, inh: str) -> str:
    """File name of the protein(-ligand) complex PDB written by senda-complex."""
    return f"{mut}_{inh}_complex.pdb"


def legacy_complex_pdb_name(mut: str, inh: str) -> str:
    """Name used before 2026-10-29 (every system was called a dimer)."""
    return f"{mut}_{inh}_dimer.pdb"


def find_complex_pdb(dirs, mut: str, inh: str) -> "Path":
    """
    First existing complex PDB in *dirs* (searched in order): the current
    {mut}_{inh}_complex.pdb, else the legacy {mut}_{inh}_dimer.pdb, so
    projects set up before the rename keep working.
    """
    from pathlib import Path

    dirs = [Path(d) for d in dirs]
    for name in (complex_pdb_name(mut, inh), legacy_complex_pdb_name(mut, inh)):
        for d in dirs:
            if (d / name).exists():
                return d / name
    raise FileNotFoundError(
        f"Complex PDB {complex_pdb_name(mut, inh)} (or legacy "
        f"{legacy_complex_pdb_name(mut, inh)}) not found in: "
        + ", ".join(str(d) for d in dirs)
    )


def protein_chains(pdb_path) -> list:
    """Chain IDs with amino-acid residues, in file order (blank chain -> 'A')."""
    chains = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[17:20].strip() not in AMINO_ACIDS:
                continue
            ch = line[21].strip() or "A"
            if ch not in chains:
                chains.append(ch)
    return chains


def resolve_chains(cfg: dict, pdb_path=None) -> list:
    """
    Chains to work on: michaelis_complex.chains if set, else every protein
    chain in *pdb_path* (monomer -> [A], tetramer -> [A, B, C, D]), else [A].
    The first entry is the default chain (no suffix on its output paths).
    """
    chains = (cfg.get("michaelis_complex") or {}).get("chains")
    if chains:
        return [str(c) for c in chains]
    if pdb_path is not None:
        found = protein_chains(pdb_path)
        if found:
            return found
    return ["A"]


def system_chains(cfg: dict, cwd, inh: str, mut: str) -> list:
    """resolve_chains() for one simulated system, reading its complex PDB if needed."""
    from pathlib import Path

    if (cfg.get("michaelis_complex") or {}).get("chains"):
        return resolve_chains(cfg)
    cwd = Path(cwd)
    protein_dir = cwd / (cfg.get("michaelis_complex") or {}).get("output_dir", "protein")
    try:
        pdb = find_complex_pdb([cwd / "simulations" / inh / mut / "replica_1" / "00_prep", protein_dir], mut, inh)
    except FileNotFoundError:
        pdb = None
    return resolve_chains(cfg, pdb)

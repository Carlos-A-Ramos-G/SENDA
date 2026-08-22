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

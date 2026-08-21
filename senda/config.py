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

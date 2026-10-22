#!/bin/bash
# Re-evaluate frames of the AMBER production with sander and with OpenMM and
# compare the energies term by term. CPU only (no GPU, no pmemd licence);
# a few minutes for ~100k atoms. Extra arguments are passed on, e.g.
#   bash 2_energy_check.sh --nframes 100 --skip-ns 5
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need python; need sander

shopt -s nullglob
trajs=("${AMBER_REPLICA}"/04_NVT/structure_NVT_[1-9]*.nc)
[ ${#trajs[@]} -gt 0 ] || { echo "ERROR: no ${AMBER_REPLICA}/04_NVT/structure_NVT_*.nc"; exit 1; }

$COMPARE energy-check --prep "${PREP}" --traj "${trajs[@]}" --out "${ENERGY_CHECK}" "$@"

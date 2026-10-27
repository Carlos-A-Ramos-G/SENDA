#!/bin/bash
# Run OpenMM stages 01_min..04_NVT on one GPU. Either run it directly
#   bash 3_run_openmm.sh
# or submit it: edit the #SBATCH lines below and `sbatch 3_run_openmm.sh`.
# Safe to resubmit: finished stages are skipped and MD resumes from the last
# checkpoint.
#SBATCH --job-name=openmm_check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gres=gpu:1
#SBATCH --time=24:00:00
#SBATCH --output=openmm_check-%j.out
##SBATCH --partition=gpu
##SBATCH --account=XXXXXXXX

# sbatch runs a spooled copy of this script, so locate the folder via the submit dir
if [ -n "${SLURM_SUBMIT_DIR:-}" ]; then cd "${SLURM_SUBMIT_DIR}"; else cd "$(dirname "${BASH_SOURCE[0]}")"; fi
[ -f ./_common.sh ] || { echo "ERROR: submit from inside example/openmm_check"; exit 1; }
source ./_common.sh
[ -x "${OPENMM_RUN}/run_openmm.sh" ] || { echo "ERROR: run 1_prepare.sh first"; exit 1; }
bash "${OPENMM_RUN}/run_openmm.sh"

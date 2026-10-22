# Settings shared by the numbered scripts in this folder. Edit, then run
# 1_prepare.sh -> 2_energy_check.sh -> 3_run_openmm.sh -> 4_compare.sh.
# Relative paths are relative to this folder.

# --- Inputs from the AMBER run you are comparing against ---------------------
# 00_prep of the AMBER replica (structure.parm7, structure_HMR.parm7,
# structure.rst7, the *_dimer.pdb). Copy it here or point to it.
PREP=./00_prep

# The AMBER replica directory itself (00_prep, 01_min ... 04_NVT, restr).
# Needed by 2_energy_check.sh (its trajectories) and 4_compare.sh (all its
# outputs); 1_prepare.sh only reads restr from it. Leave as is if you copy
# the whole replica here as ./amber_replica.
AMBER_REPLICA=./amber_replica

# NMR restraints the AMBER run used ("restr" in the replica root). Leave empty
# to take $AMBER_REPLICA/restr; set to "" and remove that file if the run had none.
RESTR=

# The config.yaml the AMBER run was set up with: the OpenMM run copies its
# protocol (HMR, temperature, restraint schedule, heating/equilibration
# lengths, output intervals) so both engines run the same thing.
CONFIG=./amber_config.yaml

# --- OpenMM run ---------------------------------------------------------------
PROD_NS=20          # production to simulate with OpenMM (and to compare)
CHUNK_NS=10         # production chunk length
PLATFORM=CUDA       # CUDA | OpenCL | CPU
SEED=2026

# Command that runs senda-omm in an environment with OpenMM + senda installed.
OMM_COMMAND="senda-omm"

# Shell lines that make python (with senda, openmm, scipy, matplotlib, pyyaml)
# and AmberTools (cpptraj, sander) available, e.g.
#   source ~/miniforge3/etc/profile.d/conda.sh && conda activate senda-omm
ENV_SETUP=""

# --- Outputs (inside this folder) ---------------------------------------------
OPENMM_RUN=./openmm_run
ENERGY_CHECK=./energy_check
REPORT=./report

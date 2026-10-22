# OpenMM check against an existing AMBER run

This folder validates the OpenMM engine (`amber_simulator.engine: openmm`) against a replica you already simulated with pmemd.cuda. pmemd.cuda is never run here. OpenMM starts from the exact same `00_prep` as the AMBER run (same parm7, coordinates and NMR restraints) and runs the same protocol. You get two independent checks:

| Step | What it tells you | Needs |
|---|---|---|
| `2_energy_check.sh` | Force-field parity: frames of the AMBER trajectory re-evaluated by `sander` and OpenMM, term by term | CPU, minutes |
| `3_run_openmm.sh` + `4_compare.sh` | Dynamics parity: heating, density, temperature, RMSD/Rg, RMSF, NMR restraints and speed vs. the AMBER run | 1 GPU, hours |

## 1. Environment (on the machine that runs the check)

```bash
conda create -n senda-omm -c conda-forge python=3.11 openmm ambertools scipy matplotlib pyyaml numpy
conda activate senda-omm
pip install -e /path/to/SENDA            # this repository, branch openMM-implementation
python -m openmm.testInstallation         # CUDA must be listed and agree with the other platforms
which cpptraj sander senda-omm            # all three must resolve
```

AmberTools from conda-forge is free and provides `cpptraj` and `sander`. pmemd.cuda is not needed.

## 2. Put the AMBER run's files here

```
example/openmm_check/
├── 00_prep/            <- copy of <AMBER replica>/00_prep
├── amber_replica/      <- the AMBER replica (00_prep, 01_min ... 04_NVT, restr), copied or symlinked
└── amber_config.yaml   <- the config.yaml that AMBER run was set up with
```

- `00_prep/` alone is enough to build and run OpenMM (steps 1 and 3).
- `2_energy_check.sh` needs `amber_replica/04_NVT/*.nc`. `4_compare.sh` needs all of `amber_replica/`'s outputs (`.out` and `.nc` of every stage).
- The NMR restraints are read from `amber_replica/restr`. OpenMM uses that exact file, so both engines are restrained identically. If you only copy `00_prep`, also copy `restr` and set `RESTR=./restr` in `settings.sh`.
- A symlink is fine for the large trajectories: `ln -s /path/to/replica_1 amber_replica`.

Everything you copy here, and everything generated here, is git-ignored.

## 3. Edit `settings.sh`

- **Paths** (`PREP`, `AMBER_REPLICA`, `RESTR`, `CONFIG`) default to the layout above.
- **`PROD_NS`:** how much production OpenMM runs, and the window compared against the start of the AMBER production. 20 ns is a reasonable first check.
- **`ENV_SETUP`:** the lines that activate the environment from step 1, e.g. `source ~/miniforge3/etc/profile.d/conda.sh && conda activate senda-omm`.
- **`PLATFORM`:** `CUDA` on the GPU node.

## 4. Run

```bash
cd example/openmm_check
bash 1_prepare.sh          # builds ./openmm_run from 00_prep (no tleap), prints the protocol
bash 2_energy_check.sh     # ./energy_check/energy_check.md   (CPU, can run right away)
sbatch 3_run_openmm.sh     # or: bash 3_run_openmm.sh  -- edit the #SBATCH lines first
bash 4_compare.sh          # ./report/summary.md + ./report/comparison.png
```

- `1_prepare.sh` prints the system, the NMR restraint source, the number of positional-restraint atoms, the equilibration cycles, the production length and the platform. Check them before submitting.
- `3_run_openmm.sh` is restartable: completed stages are skipped and MD resumes from its last checkpoint, so resubmit it if it hits the wall time.
- Extra arguments are passed through, e.g. `bash 2_energy_check.sh --nframes 100 --skip-ns 10`.

## 5. What "passing" looks like

**`energy_check.md`**
- **Bond / angle / dihedral:** |Δ| around 1e-3 kcal/mol or less.
- **Raw nonbonded:** differs by about 3.5e-5 × the electrostatic energy. That is the Coulomb-constant convention (OpenMM 332.0637 vs AMBER 18.2223² = 332.0522); the report prints the expected value.
- **Nonbonded with the constant matched:** well below 1 kcal/mol, with r ≈ 1 across frames. A larger difference means the force fields disagree, and that should be investigated before trusting any dynamics.

**`summary.md` / `comparison.png`.** These compare two independent trajectories, so expect agreement within normal replica-to-replica variation, not identical numbers:
- **Topologies:** must be identical.
- **Heating:** temperature reaches T and holds. AMBER heats with Berendsen, OpenMM with Monte Carlo pressure control, so the density during heating can differ.
- **Equilibration:** NPT density converges to the same value, within about 0.01 g/mL.
- **Production:** T, CA RMSD, Rg and the chain A–B distance in the same range, with RMSF correlated across residues.
- **NMR restraints:** inside their r2–r3 window in both engines.
- **Minimisation:** OpenMM's L-BFGS often stops before AMBER's `convergence_threshold` in mixed precision. It is then reported as "not converged" and the last cycle is used, which is expected. With the default `max_cycles_cap: 100` that can take long; lower `amber_simulator.min.max_cycles_cap` in `amber_config.yaml` (e.g. 5) for the check.

# senda

Modular Python framework for enzyme QM/MM free energy calculations.

## Installation

```bash
pip install -e .
```

**Requirements:** AmberTools >= 22 (`antechamber`, `parmchk2`, `tleap`, `pmemd.cuda` in `$PATH`).
RDKit is recommended for automatic net charge detection.

```bash
conda install -c conda-forge rdkit
```

---

## Workflow overview

```
raw crystal PDBs  ->  senda-complex  ->  curated protein-ligand dimers
                                                  |
ligand PDBs       ->  senda-param   ->  GAFF parameters
                                                  |
                       senda-sim    ->  AMBER MD replicas
                                                  |
                       senda-pmf    ->  free energy profiles  (planned)
```

On an HPC cluster, use `senda-slurm` to generate and chain all SLURM scripts automatically.

---

## Commands

| Command | What it does |
|---|---|
| `senda-param` | Parameterize ligands with GAFF2 (AM1-BCC or RESP charges) |
| `senda-complex` | Prepare curated protein-ligand dimer PDBs from raw crystal structures |
| `senda-sim` | Generate AMBER MD replica directories and optionally submit jobs |
| `senda-slurm` | Generate a chained SLURM workflow script for the full pipeline |

---

## Ligand parameterization -- `senda-param`

Parameterizes small-molecule ligands using the AMBER/GAFF force field and writes topology files for downstream MD.

### Inputs

Place one PDB file per ligand (with all hydrogens, correct geometry) in a `ligands/` folder:

```
your_project/
+-- config.yaml
+-- ligands/
    +-- LIG1.pdb
    +-- LIG2.pdb
```

### Run

```bash
senda-param --config config.yaml
```

### Which ligands are processed

Priority order:

1. `--ligands` CLI flag
2. `inhibitors:` list in `config.yaml`
3. All `.pdb` files found in `ligands/` (fallback)

### Outputs

```
ligands_libraries/
+-- LIG1/
|   +-- LIG1.mol2     # GAFF atom types + partial charges
|   +-- LIG1.frcmod   # missing GAFF parameters
|   +-- LIG1.lib      # AMBER library file
+-- LIG2/
    +-- ...
```

### Options

```
--ligands LIG1 LIG2   override the ligand list
--jobs 4              run up to 4 ligands in parallel
--force               re-run even if outputs already exist
```

---

## Michaelis complex preparation -- `senda-complex`

Converts raw crystallographic PDB files into curated protein-ligand dimer PDBs ready for AMBER MD.

### What it does

For each structure listed in `michaelis_complex.structures`:

1. **Clean** -- strips CONECT/END records and unwanted HETATM residues; keeps only the specified chains; applies `residue_renames`; copies protonation states (HID/HIE/HIP, CYM, ASH, GLH, etc.) from the reference enzyme PDB, skipping any positions that are mutated.
2. **Align** -- superimposes the cleaned structure onto the alignment reference by Ca RMSD (Kabsch). ANISOU records are carried through.
3. **Write output PDBs** -- one per inhibitor listed in `inhibitor_sources`, plus an APO (protein-only) file:
   - `{mutant}_{inhibitor}_dimer.pdb` -- protein + ligand
   - `{mutant}_APO_dimer.pdb` -- protein only

Ligand coordinates come from either the raw PDB itself (`native`) or a reference PDB you supply.

### Inputs

```
your_project/
+-- config.yaml
+-- raw_pdbs/
|   +-- wildtype.pdb
|   +-- mutant1.pdb
+-- protein/
    +-- WT_LER_dimer.pdb      # alignment reference (must exist)
    +-- WT_NIR_dimer.pdb      # NIR coordinate source (if NIR is not native)
```

### Run

```bash
senda-complex --config config.yaml
senda-complex --config config.yaml --force   # re-run even if outputs exist
```

### Filtering by mutants list

If a top-level `mutants:` list is present in the config, only structures whose mapped mutant name appears in that list are processed. Remove or leave `mutants:` empty to process all structures.

---

## AMBER MD setup -- `senda-sim`

Generates AMBER MD replica directory trees and optionally submits or launches jobs.

### Directory layout

Each replica is set up at `simulations/{inhibitor}/{mutant}/replica_{n}/`:

```
replica_1/
+-- 00_prep/   tleap inputs, topology, HMR
+-- 01_min/    minimization input
+-- 02_heat/   heating input
+-- 03_equil/  NPT + NVT equilibration inputs
+-- 04_NVT/    production inputs
+-- run_gpu    master SLURM script (cluster mode)
+-- run_local  master bash script (local mode)
```

### Run

```bash
# Generate replica directories (cluster mode, no submission)
senda-sim --config config.yaml setup

# Generate and submit immediately
senda-sim --config config.yaml setup --submit

# Local GPU workstation (no SLURM)
senda-sim --config config.yaml setup --mode local

# Submit already-generated scripts
senda-sim --config config.yaml submit
```

### Restraints

Two types of restraints are supported, both optional:

**Positional restraints** -- backbone atoms held during heating and NPT equilibration. The mask, heating weight, and per-cycle schedule are configurable.

**NMR restraints** -- distance, angle, or dihedral restraints written to an AMBER DISANG file at job time (after tleap builds the topology). Defined as a list of atom specs in the config.

---

## SLURM workflow -- `senda-slurm`

Generates a set of SLURM scripts and a top-level `submit.sh` that chains the full pipeline with `--dependency=afterok` so downstream jobs are cancelled on upstream failure.

### Job chain

```
senda-param    (optional, skip with --skip-param)
    -> afterok
senda-complex
    -> afterok
senda-launch   (runs senda-sim setup, then sbatch each replica's run_gpu)
    ->
[run_gpu x N replicas -- all run in parallel]
```

### Run

```bash
senda-slurm --config config.yaml
senda-slurm --config config.yaml --skip-param    # skip if parameters already exist
senda-slurm --config config.yaml --out-dir workflow/
```

Then submit the whole pipeline with:

```bash
bash workflow/submit.sh
```

---

## Configuration reference

All settings live in a single `config.yaml`.

```yaml
# ---- System ------------------------------------------------------------------
inhibitors:
  - LER
  - NIR
  - APO

mutants:
  - WT
  - E166V

replicas: 10

# ---- Ligand parameterization -------------------------------------------------
ligand_parameters:
  charge_method: bcc    # bcc (AM1-BCC, fast) or resp (Gaussian, accurate)
  net_charge: auto      # auto = detect via RDKit, or an explicit integer
  multiplicity: 1
  atom_type: gaff2      # gaff2 | gaff | amber

  forcefield:
    ligand: leaprc.gaff2

  # Only used when charge_method: resp
  gaussian_opt:
    nproc: 16
    mem: "32GB"
    route: "#P b3lyp/6-31g* opt"

  gaussian_hf:
    nproc: 16
    mem: "32GB"
    route: "#p hf/6-31g(d) SCF=Tight Pop=MK IOp(6/33=2)"

# ---- Michaelis complex -------------------------------------------------------
michaelis_complex:
  raw_pdbs_dir: raw_pdbs/
  output_dir: protein/
  chains: [A, B]

  # Rename residues before any other processing (e.g. numeric ligand codes).
  residue_renames:
    - from: "216"
      to: LER

  # PDB used for Ca superposition and missing-residue filling.
  alignment_reference: protein/WT_LER_dimer.pdb

  # PDB from which protonation states are read.
  # Can be holo -- ligands and water are ignored automatically.
  # Defaults to alignment_reference if omitted.
  reference_enzyme_pdb: protein/WT_LER_dimer.pdb

  # Ligand coordinate sources.
  #   native          -> already present in the raw PDB (after residue_renames)
  #   path/to/ref.pdb -> coordinates copied from this reference PDB
  inhibitor_sources:
    LER: native
    NIR: protein/WT_NIR_dimer.pdb

  # Map raw PDB filename -> mutant name used in output file names.
  structures:
    wildtype.pdb:  WT
    E166V.pdb:     E166V
    # add more mutants here

# ---- AMBER simulation --------------------------------------------------------
amber_simulator:
  use_hmr: true
  temperature: 300.0   # K

  production:
    ns_per_chunk: 100
    total_chunks: 1
    chunks_per_job: 1

  min:
    maxcyc: 30000
    ncyc: 500
    cut: 10.0
    max_cycles_cap: 100
    convergence_threshold: 3.0e-3

  heat:
    ps: 200

  equil:
    npt_ns: 1.25
    nvt_ns: 5.0

  output:
    ntpr: 10000
    ntwx: 500000
    ntwr: 10000

  restraints:
    # Positional (backbone) restraints during heating and NPT equilibration.
    positional:
      mask: "@CA,C,N,O,H &!:WAT"
      heating_weight: 20.0                           # kcal/mol/A^2
      equil_schedule: [15.0, 12.0, 9.0, 6.0, 3.0]  # one value per NPT cycle

    # NMR restraints (distance, angle, dihedral). Remove if not needed.
    # atoms: list of {residue: <resname>, name: <atomname>}
    # 2 atoms = distance, 3 = angle, 4 = dihedral
    # Add index: <n> (0-based) for cross-residue atoms with multiple matches.
    nmr:
      - type: dihedral
        atoms:
          - {residue: LER, name: O5}
          - {residue: LER, name: C7}
          - {residue: LER, name: C6}
          - {residue: LER, name: N3}
        r1: -10.0
        r2:  -5.0
        r3:  15.0
        r4:  20.0
        rk2: 500
        rk3: 500

# ---- SLURM -------------------------------------------------------------------
slurm:
  amber_module: apps/amber/24

  # Command to activate the Python environment where senda is installed.
  # Examples: "conda activate senda" | "module load python/3.11"
  senda_env: ""

  param:          # ligand parameterization (CPU)
    time: "4:00:00"
    ntasks: 1
    cpus-per-task: 4
    mem: "8G"
    partition: cpu
    account: MY_ACCOUNT

  complex:        # senda-complex (CPU)
    time: "2:00:00"
    ntasks: 1
    cpus-per-task: 1
    mem: "8G"
    partition: cpu
    account: MY_ACCOUNT

  launcher:       # senda-sim setup + sbatch all run_gpu scripts (CPU)
    time: "0:30:00"
    ntasks: 1
    cpus-per-task: 1
    mem: "4G"
    partition: cpu
    account: MY_ACCOUNT

  master:         # stages 00-03 (GPU)
    time: "1-00:00:00"
    ntasks: 1
    gres: gpu:1
    partition: gpu
    account: MY_ACCOUNT

  nvt:            # NVT production chunks (GPU)
    time: "5-00:00:00"
    ntasks: 1
    gres: gpu:1
    partition: gpu
    account: MY_ACCOUNT
```

---

## Charge methods

| Method | Speed | Accuracy | Requires |
|---|---|---|---|
| `bcc` | Fast | Good | AmberTools only |
| `resp` | Slow | High | AmberTools + Gaussian on HPC |

**BCC** runs entirely locally: `antechamber -> parmchk2 -> tleap`.

**RESP** is a two-phase process. Phase 1 generates Gaussian input files locally. Submit those to HPC, then re-run `senda-param` to complete the parameterization. The command prints the exact steps at each phase.

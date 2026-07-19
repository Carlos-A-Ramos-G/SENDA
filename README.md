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
                       senda-analyse ->  representative frame (rst7)
                                                  |
                       senda-qmmm string equil   ->  QM/MM equilibration
                                                  |
                       senda-qmmm string prod    ->  restraint-free production (optional)
                                                  |
                       senda-qmmm string scan    ->  restrained path scan
                                                  |
                       senda-qmmm string string  ->  free energy profile (ASM)
```

On an HPC cluster, use `senda-slurm` to generate and chain all SLURM scripts automatically.

---

## Commands

| Command | What it does |
|---|---|
| `senda-param` | Parameterize ligands with GAFF2 (AM1-BCC or RESP charges) |
| `senda-complex` | Prepare curated protein-ligand dimer PDBs from raw crystal structures |
| `senda-sim` | Generate AMBER MD replica directories, submit jobs, check run status, or inspect topology |
| `senda-slurm` | Generate a chained SLURM workflow script for the full pipeline |
| `senda-analyse` | Analyse reactive distances from NVT trajectories and select a representative frame for QM/MM |
| `senda-qmmm string equil` | Stage 05: set up and optionally submit the QM/MM equilibration job |
| `senda-qmmm string prod` | Stage 05_QMMM_restraint_free: optional unrestrained QM/MM production run after equilibration |
| `senda-qmmm string scan` | Stage 06: set up the restrained scan along the initial guess path |
| `senda-qmmm string string` | Stage 07: set up the adaptive string method (ASM) calculation |

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

### Config blocks read

`senda-complex` reads three blocks from `config.yaml`:

| Block | Role |
|---|---|
| `michaelis_complex:` | All structural settings (required) |
| `inhibitors:` | Top-level list used to filter which inhibitors are written (optional — all if absent) |
| `mutants:` | Top-level list used to filter which structures are processed (optional — all if absent) |

Everything else in the config (`amber_simulator:`, `analysis:`, `qmmm:`, etc.) is ignored.

### What it does

For each structure in `michaelis_complex.structures` whose mutant name appears in the top-level `mutants:` list (or all structures if `mutants:` is absent or empty):

1. **Clean** -- strips CONECT/END records and unwanted HETATM residues; keeps only the specified chains; applies `residue_renames`; copies protonation states (HID/HIE/HIP, CYM, ASH, GLH, etc.) from the reference enzyme PDB, skipping any positions that are mutated.
2. **Align** -- superimposes the cleaned structure onto the alignment reference by Ca RMSD (Kabsch). ANISOU records are carried through.
3. **Fill missing residues** -- residues present in the alignment reference but absent in the raw PDB are copied from the reference (valid because the structure has already been superimposed).
4. **Write output PDBs** -- one per inhibitor listed in `inhibitor_sources`, plus an APO (protein-only) file:
   - `{mutant}_{inhibitor}_dimer.pdb` -- protein + ligand
   - `{mutant}_APO_dimer.pdb` -- protein only

Ligand coordinates come from either the raw PDB itself (`native`) or a reference PDB you supply.

### Inputs

The alignment reference PDB must already exist before running `senda-complex`. Typically this is a manually curated WT + ligand structure that serves as the common coordinate frame for all other structures.

```
your_project/
+-- config.yaml
+-- raw_pdbs/
|   +-- wildtype.pdb
|   +-- mutant1.pdb
+-- protein/
    +-- WT_LER_dimer.pdb      # alignment reference -- must exist beforehand
    +-- WT_NIR_dimer.pdb      # NIR coordinate source (if NIR is not native)
```

### Run

```bash
senda-complex --config config.yaml
senda-complex --config config.yaml --force   # re-run even if outputs exist
```

### Filtering by inhibitors and mutants

`senda-complex` respects both top-level lists:

- **`inhibitors:`** — only inhibitors listed here are prepared. `APO` is a special value: include it to also write protein-only (APO) PDB files; omit it to skip APO output entirely. If `inhibitors:` is absent or empty, all entries in `michaelis_complex.inhibitor_sources` are processed and APO files are always written.
- **`mutants:`** — only structures whose mapped mutant name appears here are processed. Remove or leave empty to process all structures.

### Residue trimming

After alignment, each chain is clipped to the residue number range of the alignment reference. This prevents crystal structures that resolve extra C-terminal (or N-terminal) residues from introducing spurious charge differences between systems.

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

### Check simulation status

```bash
# Report progress for every replica
senda-sim --config config.yaml check

# Show only replicas that are not fully complete
senda-sim --config config.yaml check --failed

# Also print the last lines of the in-progress output file
senda-sim --config config.yaml check --failed -v
```

Output shows two columns per replica: **Completed** (last stage with the AMBER completion marker) and **In progress** (stage currently running or interrupted). A truncated `.out` file is never labelled "failed" — it may belong to a still-running job.

```
  System               Completed              In progress
  -------------------------------------------------------
  LER/WT/replica_1     03_equil  6/6          04_NVT  3/20
  LER/WT/replica_2     04_NVT   20/20         --
  LER/WT/replica_3     02_heat  200000 steps  03_equil  2/6
  LER/WT/replica_4     --                     --
```

### Topology info

```bash
senda-sim --config config.yaml topology_info
```

Reads `replica_1/00_prep/structure.parm7` for each (inhibitor, mutant) pair and prints a table of atom count, box dimensions, water count, Na+, Cl-, and the net charge of the solute (protein + ligand, excluding water and ions). Useful for verifying that ion placement is consistent across systems. Requires ParmEd (`pip install parmed` or load AmberTools).

### Simulation protocol

The MD protocol applied to each replica is:

| Stage | Type | Details |
|---|---|---|
| `00_prep` | tleap | Solvation, ion placement, HMR |
| `01_min` | Minimization | Convergence-based loop, up to `max_cycles_cap` |
| `02_heat` | NPT (Berendsen) | 1 K to target temperature; retry loop on crash |
| `03_equil` | NPT (MC) -> NVT | 5 restrained NPT cycles + 1 unrestrained NVT |
| `04_NVT` | NVT | Production, chained chunks |

Heating uses the Berendsen barostat (`barostat=1`) for stability during the far-from-equilibrium temperature ramp, then switches to the Monte Carlo barostat for equilibration. If heating crashes, the run script retries automatically (up to `heat.max_retries` attempts), restarting from the last written checkpoint if one exists.

### Ion placement -- SPLIT method

By default, senda neutralises the system with `addions Na+ 0`. To add NaCl at a specific concentration using the SPLIT method (Machado & Pantano, *J. Chem. Theory Comput.* 2020), set `leap.salt_molarity` in the config:

```yaml
amber_simulator:
  leap:
    salt_molarity: 0.150   # mol/L
```

This triggers a two-pass tleap scheme: pass 1 determines the solute charge and the number of solvation waters; the run script computes the exact Na+/Cl- counts at the target molarity and substitutes them before pass 2 builds the final system.

### Restraints

Two types of restraints are supported, both optional:

**Positional restraints** -- backbone atoms held during heating and NPT equilibration. The mask, heating weight, and per-cycle schedule are configurable.

**NMR restraints** -- distance, angle, or dihedral restraints written to an AMBER DISANG file at job time (after tleap builds the topology). Must be a dict keyed by inhibitor name; each inhibitor's restraints are applied only to that inhibitor's topology. Omit an inhibitor's key (or set it to `[]`) for no NMR restraints for that ligand.

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

## NVT trajectory analysis -- `senda-analyse`

Analyses reactive distances from pooled NVT trajectories across all replicas and selects a single representative frame as the starting point for QM/MM relaxation.

### What it does

For each inhibitor/mutant combination:

1. Loads all `04_NVT/structure_NVT_*.nc` frames from every replica into a pooled dataset
2. Computes the probability distribution of each reactive distance
3. Identifies the most probable value (mode) for each distance; warns if the distribution is multimodal
4. If `water_rdf` is configured: computes the radial distribution function g(r) of WAT-O atoms around the specified center, finds the first coordination shell peak, and adds a soft score penalty to frames that have no water within peak +/- tolerance
5. Scores every frame by the sum of sigma-normalised distance deviations plus any water penalties
6. Selects the frame where all distances fall within mode +/- 1 sigma (relaxes to 2 or 3 sigma if needed), preferring frames that also satisfy the water criterion
7. Writes the selected frame as a restart file, distance distribution plots, and (if configured) RDF plots

### Outputs

One set of files per chain (A, B, ... from `michaelis_complex.chains`):

```
simulations/{inhibitor}/{mutant}/
+-- {inhibitor}_{mutant}_chainA_representative.rst7
+-- {inhibitor}_{mutant}_chainB_representative.rst7

Analysis/
+-- data/
|   +-- distances/
|   |   +-- {inhibitor}_{mutant}_chainA_distances.csv
|   |   +-- {inhibitor}_{mutant}_chainB_distances.csv
|   +-- rdf/                              (only if water_rdf is configured)
|       +-- {inhibitor}_{mutant}_chainA_{label}.csv
|       +-- {inhibitor}_{mutant}_chainB_{label}.csv
+-- plots/
    +-- distances/
    |   +-- {inhibitor}_{mutant}_chainA_distances.png
    |   +-- {inhibitor}_{mutant}_chainB_distances.png
    +-- rdf/                              (only if water_rdf is configured)
        +-- {inhibitor}_{mutant}_chainA_{label}.png
        +-- {inhibitor}_{mutant}_chainB_{label}.png
```

Distance CSV columns: `replica,frame,<label1>,<label2>,...` (frame is 1-based within each replica's concatenated trajectory; VMD frame = CSV frame - 1).

RDF CSV columns: `r,g_r`.

### Run

```bash
senda-analyse --config config.yaml
senda-analyse --config config.yaml --inhibitors LER NIR   # subset of inhibitors
senda-analyse --config config.yaml --mutants H172Y        # subset of mutants
```

APO systems are skipped automatically (reactive distances require a ligand).

### Specifying atoms

Each atom in `reactive_distances` and `water_rdf.center_atoms` uses one of two fields to identify its residue:

| Field | When to use | Value |
|---|---|---|
| `sequence` | Protein / binding-site residues | PDB residue number as it appears in the dimer PDB (same in every monomer) |
| `substrate_sequence` | Ligand or substrate residues | 1-based AMBER residue number of the substrate in `chains[0]`, as assigned in the parm7 topology |

For `sequence` atoms the equivalent residue in every other monomer is looked up directly from `chain_map[(chain, pdb_resnum)]`. This works as long as protein chains share PDB residue numbers, which is standard for homodimers and homo-oligomers.

For `substrate_sequence` atoms the equivalent residue in every other monomer is determined by **geometric proximity**: the candidate non-solvent residue whose atoms are closest to the declared `sequence` reference atoms (resolved to that chain) is selected. This is necessary because ligands may share a PDB residue number with a protein residue (making direct lookup ambiguous) and because a peptide substrate uses standard amino-acid names that cannot be distinguished by name alone.

To find the correct AMBER resid for your ligand/substrate in `chains[0]`, inspect the parm7 topology with `parmed` or `pytraj`.

### Water RDF

`water_rdf` is an optional list of radial distribution function specifications. Each entry computes g(r) of WAT-O atoms around a reference point and uses the result to add a soft score penalty to frames that have no water near the first coordination shell peak.

```yaml
water_rdf:
  - label: "water_CYS145_SG"   # used in output file names
    center_atoms:               # one or more atoms; multiple atoms use their centroid
      - {sequence: 145, name: SG}
    r_max: 10.0                 # maximum radius in Angstroms (default 10.0)
    tolerance: 0.3              # half-width of the peak window in Angstroms (default 0.3)
    penalty: 3.0                # score penalty when no water in window (default 3.0)
```

The penalty is added to the distance-deviation score (lower = better), so a value of 3.0 is equivalent to one distance being 3 sigma from its mode. Frames that satisfy the distance criterion but lack a water at the first peak are deprioritised rather than excluded.

### Requirements

```bash
pip install senda[analysis]   # installs scipy, matplotlib, pytraj
```

---

## QM/MM string method -- `senda-qmmm string`

Drives the QM/MM adaptive string method (ASM) workflow using AMBER's `sander.MPI`. Each stage is independently callable so individual steps can be re-run without restarting the full pipeline.

### Prerequisites

```bash
pip install senda[qmmm]   # installs parmed, pytraj
```

AMBER with `sander.MPI` and `cpptraj` must be in `$PATH` (or loaded via a module on the cluster).

### Stages

| Subcommand | Stage | What it does |
|---|---|---|
| `senda-qmmm string equil` | 05 | Resolves CV atoms and QM region, builds H10 topology, writes AMBER input + SLURM script for QM/MM equilibration |
| `senda-qmmm string prod` | 05_QMMM_restraint_free | **Optional.** Writes AMBER input + SLURM script for an unrestrained QM/MM production run starting from the equil output. If run, the scan stage uses its restart file as the starting structure. |
| `senda-qmmm string scan` | 06 | Writes per-node harmonic restraint files from the interpolated guess, AMBER input template, scan SLURM script, and cpptraj centering script |
| `senda-qmmm string string` | 07 | Writes CVs file, string guess, per-node input files (`in.sh`), groupfile, and SLURM script for `sander.MPI -ng N -groupfile` |

### Run

```bash
# Write all input files (no job submission)
senda-qmmm string equil  config.yaml
senda-qmmm string prod   config.yaml   # optional: restraint-free production
senda-qmmm string scan   config.yaml
senda-qmmm string string config.yaml

# Write and submit to SLURM
senda-qmmm string equil  -s config.yaml
senda-qmmm string prod   -s --after <equil_jobid>  config.yaml   # optional
senda-qmmm string scan   -s --after <equil_or_prod_jobid>  config.yaml
senda-qmmm string string -s --after <scan_jobid>   config.yaml

# Process only one inhibitor / mutant
senda-qmmm string equil -i LER -m WT config.yaml
```

### Outputs

All output lives alongside the replica directories for each inhibitor/mutant combination:

```
simulations/{inhibitor}/{mutant}/
+-- structure_H10.parm7           # topology with QM hydrogen masses set to 10 amu
+-- _guess_interpolated.npy       # arc-length-interpolated guess (internal cache)
+-- _qmmm_string_meta.json        # resolved metadata shared across stages
+-- 05_QMMM_equilibration/
|   +-- in                        # AMBER QM/MM input
|   +-- restr                     # extra restraints (AMBER &rst blocks)
|   +-- equilibration.sh          # SLURM script
+-- 05_QMMM_restraint_free/       # only present if senda-qmmm string prod was run
|   +-- in                        # AMBER QM/MM input (no restraints, irest=1)
|   +-- prod.sh                   # SLURM script
+-- 06_QMMM_scan/
|   +-- in_template               # AMBER input with __NODE__ placeholder
|   +-- restr0                    # extra restraints appended per node by scan job
|   +-- restr{1..N}               # per-node CV harmonic restraints
|   +-- scan.sh                   # SLURM script (sequential node loop)
|   +-- center.sh                 # cpptraj centering script (run at end of scan)
+-- 07_QMMM_string/
    +-- in                        # AMBER string input (__SEED__ filled by in.sh)
    +-- in.sh                     # generates per-node {i}.in files + string.groupfile
    +-- guess                     # string guess with AMBER header (N_nodes  N_cvs  0.0)
    +-- CVs                       # AMBER CVs file for sander ASM
    +-- string.sh                 # SLURM script (sander.MPI -ng N -groupfile)
```

### Collective variable atom specs

CV atoms use the same residue-by-name syntax as `senda-analyse`, plus an additional spec for catalytic water molecules:

| Spec | When to use |
|---|---|
| `{sequence: N, name: atomname}` | Protein residue by PDB residue number |
| `{substrate_sequence: N, name: atomname}` | Ligand/substrate by AMBER resid in chain A |
| `{nearest_water_to: <ref_spec>, name: atomname}` | WAT molecule whose O is closest to `ref_spec` in the representative frame |

CV types supported: `distance`, `angle`, `dihedral` (map to AMBER `BOND`, `ANGLE`, `TORSION`).

### QM region

If `qmmask` is not set, senda selects the QM region automatically:

1. Seed: all residues containing a CV atom.
2. Expand via bond-graph (parmed) until every QM/MM boundary bond is a C-C bond.
3. Protein backbone C-N peptide bonds are handled specially to ensure a valid cut.
4. Net charge is estimated by summing parmed partial charges of QM atoms.

Override by setting `qmmask` and `qmcharge` explicitly in the inhibitor config block.

### H10 topology

Hydrogen atoms in the CV definitions have their mass set to 10 amu in a modified topology (`structure_H10.parm7`). This improves sampling of light-atom CVs in the ASM. For WAT residues both hydrogens are patched even if only one appears in a CV.

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

  leap:
    salt_molarity: 0.150   # NaCl concentration (mol/L); uses SPLIT method (Machado & Pantano, JCTC 2020)
                           # omit or set to 0 for neutralisation-only (addions Na+ 0)

  heat:
    ps: 200
    max_retries: 5         # retry heating on crash, restarting from last checkpoint

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

    # NMR restraints (distance, angle, dihedral). Must be a dict keyed by
    # inhibitor name -- each inhibitor's restraints are applied only to that
    # inhibitor's topology. Omit the key or set it to [] for no restraints.
    # atoms: list of {residue: <resname>, name: <atomname>}
    # 2 atoms = distance, 3 = angle, 4 = dihedral
    # Add index: <n> (0-based) for cross-residue atoms with multiple matches.
    nmr:
      LER:
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
      NIR: []   # no NMR restraints for NIR

# ---- Analysis ----------------------------------------------------------------
analysis:
  # Reactive distances are declared per inhibitor so each ligand can track
  # its own reactive atoms independently.
  #
  # sequence:           PDB residue number as it appears in the dimer PDB.
  #                     Same in every monomer; used for protein residues.
  # substrate_sequence: AMBER residue number of the substrate in chains[0],
  #                     as assigned in the parm7 topology (check with parmed/pytraj).
  #                     Equivalent residues in other monomers are found by geometric proximity.
  # The analysis runs automatically for every chain in michaelis_complex.chains.

  LER:
    reactive_distances:
      - label: "C5-OG_SER144"   # label used in plots and output file names
        atoms:
          - {substrate_sequence: 301, name: C5}  # AMBER resid of LER in chains[0] (check parm7)
          - {sequence: 144, name: OG}            # PDB resnum of SER144
      - label: "C5-NE2_HIS41"
        atoms:
          - {substrate_sequence: 301, name: C5}
          - {sequence: 41,  name: NE2}
    water_rdf:
      - label: "water_SER144_OG"
        center_atoms:
          - {sequence: 144, name: OG}
        r_max: 10.0      # Angstroms (default 10.0)
        tolerance: 0.3   # peak +/- window (default 0.3)
        penalty: 3.0     # score penalty when no water at first peak (default 3.0)

  NIR:
    reactive_distances:
      - label: "C5-OG_SER144"
        atoms:
          - {substrate_sequence: 301, name: C5}  # AMBER resid of NIR in chains[0] (check parm7)
          - {sequence: 144, name: OG}
      - label: "C5-NE2_HIS41"
        atoms:
          - {substrate_sequence: 301, name: C5}
          - {sequence: 41,  name: NE2}
    # water_rdf: []  # omit the key or leave empty for no water RDF

# ---- SLURM -------------------------------------------------------------------
slurm:
  amber_module: apps/amber/24
  account: MY_ACCOUNT

  # Command to activate the Python environment where senda is installed.
  # Examples: "conda activate senda" | "module load python/3.11"
  senda_env: ""

  cpu:            # param, complex, and launcher jobs
    partition: cpu
    ntasks: 1
    cpus-per-task: 4
    mem: "8G"
    time: "4:00:00"

  gpu:            # run_gpu (stages 00-03) and NVT production chunks
    partition: gpu
    ntasks: 1
    gres: gpu:1
    time: "5-00:00:00"

  qmmm:           # senda-qmmm string jobs (equil and scan share ntasks; string scales automatically)
    ntasks: 8               # MPI tasks for equil + scan jobs, and tasks-per-node for string
    time: "1-00:00:00"      # wall time for equil and scan
    time_string: "7-00:00:00"  # wall time for the string method job

# ---- QM/MM string method -----------------------------------------------------
qmmm:
  string:
    inhibitors:
      LER:
        mutants: [WT, E166V]   # subset of top-level mutants to run; omit key for all

        collective_variables:
          - type: distance        # distance | angle | dihedral
            atoms:
              - {sequence: 41, name: NE2}              # HIS41 NE2 (protein)
              - {substrate_sequence: 301, name: C1}    # substrate C1 (ligand)
          - type: distance
            atoms:
              - {sequence: 145, name: SG}              # CYS145 SG (protein)
              - {nearest_water_to: {sequence: 41, name: NE2}, name: O}  # catalytic water O

        guess: guesses/LER_path.dat   # initial path; interpolated to n_nodes automatically

        qm_theory: DFTB3   # semiempirical level (DFTB3 | PM6 | AM1 | etc.)
        qmcut: 12.0         # QM electrostatic cutoff in Angstroms

        # Optional: override automatic QM region selection.
        # If set, qmcharge must also be provided.
        # qmmask: "@1-50,301-310"
        # qmcharge: -1

        extra_restraints:   # optional; appended to every restraint file
          - atoms: [12, 34]
            r1: 1.0
            r2: 2.0
            r3: 2.5
            r4: 4.0
            rk2: 50.0
            rk3: 50.0

        equil:              # stage 05 -- QM/MM equilibration
          temp: 300.0
          nstlim: 20000
          dt: 0.001
          gamma_ln: 5.0
          ntpr: 50
          ntwx: 100
          ntwr: 100

        prod:               # stage 05_QMMM_restraint_free -- optional unrestrained production
          nstlim: 100000    # run length (default 100000 steps = 100 ps at dt=0.001)
          dt: 0.001
          gamma_ln: 1.0     # lower friction than equil; typical for production
          ntpr: 500
          ntwx: 500
          ntwr: 500
          # temp: 300.0     # inherits from equil.temp if omitted

        scan:               # stage 06 -- restrained scan
          n_nodes: 64             # number of windows; must equal string.n_nodes
          force_constant: 100.0   # harmonic force constant (kcal/mol/A^2 or /rad^2)
          nstlim: 5000
          dt: 0.001
          gamma_ln: 5.0

        string:             # stage 07 -- adaptive string method
          n_nodes: 64
          nstlim: 50000
          dt: 0.001
          gamma_ln: 5.0
          seed: 1234              # base random seed; each node gets seed + node_index
          prep_steps: 500         # ASM preparation steps before string update
          z_bias: false           # Fortran logical (.false. / .true.)
          force_constant_d: 100.0 # string force constant
```

---

## Charge methods

| Method | Speed | Accuracy | Requires |
|---|---|---|---|
| `bcc` | Fast | Good | AmberTools only |
| `resp` | Slow | High | AmberTools + Gaussian on HPC |

**BCC** runs entirely locally: `antechamber -> parmchk2 -> tleap`.

**RESP** is a two-phase process. Phase 1 generates Gaussian input files locally. Submit those to HPC, then re-run `senda-param` to complete the parameterization. The command prints the exact steps at each phase.

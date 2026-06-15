# senda

Modular Python framework for enzyme QM/MM free energy calculations.

## Installation

```bash
pip install -e .
```

**Requirements:** AmberTools ≥ 22 (`antechamber`, `parmchk2`, `tleap` in `$PATH`). RDKit is recommended for automatic net charge detection.

```bash
conda install -c conda-forge rdkit biopython
```

---

## Workflow overview

```
raw crystal PDBs  →  senda-complex  →  curated dimers
                                            ↓
ligand PDBs       →  senda-param   →  GAFF parameters
                                            ↓
                      senda-sim    →  AMBER MD replicas   (planned)
                                            ↓
                      senda-fep    →  free energy profiles (planned)
```

---

## Ligand parameterization — `senda-param`

`senda-param` parameterizes small-molecule ligands using the AMBER/GAFF force field and writes the topology files needed for downstream MD simulations.

### Inputs

Place one PDB file per ligand (with all hydrogens, correct geometry) in a `ligands/` folder:

```
your_project/
├── config.yaml
└── ligands/
    ├── LIG1.pdb
    └── LIG2.pdb
```

### Run

```bash
cd your_project/
senda-param --config config.yaml
```

### Which ligands are processed

`senda-param` resolves the ligand list in this order of priority:

1. `--ligands` CLI flag — explicit override, always takes precedence
2. `inhibitors:` key in `config.yaml` — used when present
3. All `.pdb` files found in `ligands/` — fallback when neither is set

If the list comes from `inhibitors:` or `--ligands`, every entry is validated against the `ligands/` folder before any work starts. A missing PDB causes an immediate error.

### Outputs

Three files per ligand are written to `ligands_libraries/<ligand>/`:

```
ligands_libraries/
├── LIG1/
│   ├── LIG1.mol2     # GAFF atom types + partial charges
│   ├── LIG1.frcmod   # missing GAFF parameters
│   └── LIG1.lib      # AMBER library file
└── LIG2/
    ├── ...
```

### Options

```
senda-param --config config.yaml --ligands LIG1 LIG2   # override ligand list
senda-param --config config.yaml --jobs 4              # run up to 4 ligands in parallel
senda-param --config config.yaml --force               # re-run even if outputs already exist
```

---

## Michaelis complex preparation — `senda-complex`

`senda-complex` converts raw crystallographic PDB files into curated protein–ligand dimer PDBs ready for AMBER MD setup.

### What it does

For each raw crystal structure listed in the config:

1. **Clean** — removes CONECT records, glycerol (GOL), non-A/B chains (e.g. tetramer C/D chains); renames HIS residues to HID/HIE based on protonation state; renames the inhibitor residue `216` → `LER`.
2. **Align** — superimposes the cleaned structure onto a reference frame (WT_LER) by Cα RMSD minimisation (Biopython/Kabsch). ANISOU thermal ellipsoids are carried through unchanged.
3. **Write three output PDBs** per mutant:
   - `{mutant}_LER_dimer.pdb` — protein + LER from the crystal structure
   - `{mutant}_NIR_dimer.pdb` — protein + NIR placed at the canonical position (taken from `WT_NIR_dimer.pdb`, which is already in the reference frame)
   - `{mutant}_APO_dimer.pdb` — protein only

Single-chain crystals (e.g. only chain A) are handled transparently: NIR is inserted only on the chains that are present.

### Inputs

```
your_project/
├── config.yaml
├── raw_pdbs/
│   ├── wildtype_1216_refmac6.pdb
│   ├── E166V_1216_3.1_refmac17.pdb
│   └── ...
└── protein/
    ├── WT_LER_dimer.pdb     # alignment reference (must exist)
    └── WT_NIR_dimer.pdb     # NIR coordinate source (must exist)
```

### Run

```bash
cd your_project/
senda-complex --config config.yaml
```

### Outputs

```
protein/
├── WT_LER_dimer.pdb
├── WT_NIR_dimer.pdb
├── WT_APO_dimer.pdb
├── E166V_LER_dimer.pdb
├── E166V_NIR_dimer.pdb
├── E166V_APO_dimer.pdb
└── ...
```

### Options

```
senda-complex --config config.yaml           # skip already-done mutants
senda-complex --config config.yaml --force   # re-run everything
```

---

## Configuration

All workflow settings live in a single `config.yaml` in the working directory.

### Charge method (`senda-param`)

```yaml
ligand_parameters:
  charge_method: bcc    # bcc (default) or resp
```

| Method | Speed | Accuracy | Requires |
|--------|-------|----------|----------|
| `bcc` | Fast | Good | AmberTools only |
| `resp` | Slow | High | AmberTools + Gaussian on HPC |

**BCC** runs locally: `antechamber → parmchk2 → tleap`.

**RESP** generates Gaussian input files locally (Phase 1), which must be submitted to HPC. `senda-param` prints the exact commands at each step.

### Michaelis complex section (`senda-complex`)

```yaml
michaelis_complex:

  raw_pdbs_dir: raw_pdbs/       # directory with raw crystal structures
  output_dir: protein/          # where curated dimers are written

  reference_ler: protein/WT_LER_dimer.pdb   # Cα superposition target
  reference_nir: protein/WT_NIR_dimer.pdb   # canonical NIR coordinates

  # HIS protonation by residue number (both chains A and B)
  his_rename:
    41:  HID
    64:  HIE
    80:  HID
    163: HIE
    164: HIE
    172: HIE
    246: HIE

  # raw PDB filename → mutant name used in output file names
  structures:
    wildtype_1216_refmac6.pdb:   WT
    E166V_1216_3.1_refmac17.pdb: E166V
    # ... add more mutants here
```

### Full configuration reference

```yaml
inhibitors:
  - NIR
  - LER
mutants:
  - WT
  - E166V
replicas: 10

ligand_parameters:
  charge_method: bcc       # bcc | resp
  net_charge: auto         # auto (RDKit) or an explicit integer
  multiplicity: 1
  atom_type: gaff2         # gaff2 | gaff | amber

  forcefield:
    ligand: leaprc.gaff2

  # Gaussian settings — only used when charge_method: resp
  gaussian_opt:
    nproc: 16
    mem: "32GB"
    route: "#P b3lyp/6-31g* opt"

  gaussian_hf:
    nproc: 16
    mem: "32GB"
    route: "#p hf/6-31g(d) SCF=Tight Pop=MK IOp(6/33=2)"

michaelis_complex:
  raw_pdbs_dir: raw_pdbs/
  output_dir: protein/
  reference_ler: protein/WT_LER_dimer.pdb
  reference_nir: protein/WT_NIR_dimer.pdb
  his_rename:
    41: HID
    64: HIE
    # ...
  structures:
    wildtype_1216_refmac6.pdb: WT
    # ...
```

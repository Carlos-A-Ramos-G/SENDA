# senda

Modular Python framework for enzyme QM/MM free energy calculations.

## Installation

```bash
pip install -e .
```

**Requirements:** AmberTools ≥ 22 (`antechamber`, `parmchk2`, `tleap` in `$PATH`). RDKit is recommended for automatic net charge detection.

```bash
conda install -c conda-forge rdkit
```

---

## Ligand parameterization

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

These three files are what `tleap` needs to build a complex topology.

### Options

```
senda-param --config config.yaml --ligands LIG1 LIG2   # override ligand list
senda-param --config config.yaml --jobs 4              # run up to 4 ligands in parallel
senda-param --config config.yaml --force               # re-run even if outputs already exist
```

---

## Configuration

Workflow settings are controlled by a `config.yaml` file placed in the working directory. All keys are optional — sensible defaults are used when the file is absent.

### Charge method

The most important choice is how partial charges are assigned:

```yaml
ligand_parameters:
  charge_method: bcc    # bcc (default) or resp
```

| Method | Speed | Accuracy | Requires |
|--------|-------|----------|----------|
| `bcc` | Fast | Good | AmberTools only |
| `resp` | Slow | High | AmberTools + Gaussian on HPC |

**BCC** runs the full pipeline locally: `antechamber → parmchk2 → tleap`.

**RESP** generates the Gaussian input files locally (Phase 1), which must then be submitted to HPC. `senda-param` prints the exact commands to run at each step after Gaussian completes.

### Full configuration reference

```yaml
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
```

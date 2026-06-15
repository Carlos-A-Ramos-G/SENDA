"""
ligand_parameters.slurm
Generate a single SLURM batch script that runs the entire RESP parameterization
pipeline end-to-end for one ligand without manual intervention between phases:

  Phase 1  – ligand_parameters-run  (write Gaussian inputs + RESP files)
  Phase 2a – g16 < {lig}_opt.com > {lig}_opt.log   (geometry optimisation)
  Phase 2b – ligand_parameters-hf-input {lig}_opt.log        (build HF/ESP input)
  Phase 2c – g16 < {lig}_hf.com  > {lig}_hf.log    (HF single-point for ESP)
  Phase 3  – ligand_parameters-amber {lig}_hf.log             (AMBER parameterization)
"""

from pathlib import Path


_TEMPLATE = """\
#!/bin/bash
{sbatch_directives}
# -- Paths (baked in at generation time) --------------------------------------
PROJ_ROOT="{proj_root}"
LP_WD="{workdir}"
LP_CONFIG="{config_path}"
LP_LIG="{lig}"

# -- Environment ---------------------------------------------------------------
{module_lines}
export GAUSS_SCRDIR="$LP_WD"
{conda_block}
cd "$PROJ_ROOT"

# -- Phase 1: write Gaussian inputs and RESP files ----------------------------
echo "[$(date '+%H:%M:%S')] Phase 1 -- ligand_parameters-run"
ligand_parameters-run --ligand "$LP_LIG" --config "$LP_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 1" && exit 1

# -- Phase 2a: geometry optimisation ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 2a -- geometry optimisation ({opt_com})"
g16 < "$LP_WD/{opt_com}" > "$LP_WD/{opt_log}"
[ $? -ne 0 ] && echo "ERROR in Phase 2a (Gaussian opt)" && exit 1

# -- Phase 2b: build HF/ESP input from optimised geometry ---------------------
echo "[$(date '+%H:%M:%S')] Phase 2b -- ligand_parameters-hf-input"
ligand_parameters-hf-input "$LP_WD/{opt_log}" --config "$LP_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 2b (ligand_parameters-hf-input)" && exit 1

# -- Phase 2c: HF/6-31G(d) single-point for ESP/RESP -------------------------
echo "[$(date '+%H:%M:%S')] Phase 2c -- HF single-point ({hf_com})"
g16 < "$LP_WD/{hf_com}" > "$LP_WD/{hf_log}"
[ $? -ne 0 ] && echo "ERROR in Phase 2c (Gaussian HF)" && exit 1

# -- Phase 3: AMBER parameterization ------------------------------------------
echo "[$(date '+%H:%M:%S')] Phase 3 -- ligand_parameters-amber"
ligand_parameters-amber "$LP_WD/{hf_log}" --config "$LP_CONFIG"
[ $? -ne 0 ] && echo "ERROR in Phase 3 (ligand_parameters-amber)" && exit 1

echo "[$(date '+%H:%M:%S')] Done. Output files in $LP_WD"
"""


def write_slurm(lig, cfg, output, proj_root, workdir, config_path=None):
    """
    Generate the SLURM batch script for the full RESP pipeline for one ligand.

    Parameters
    ----------
    lig         : str        -- ligand ID (PDB filename stem)
    cfg         : dict       -- full parsed config.yaml
    output      : str | Path -- path of the generated script
    proj_root   : str | Path -- absolute project root (cd here before running)
    workdir     : str | Path -- absolute path to output_dir/{lig}/
    config_path : str | Path -- absolute config file path baked into the script
    """
    g_opt    = cfg.get('gaussian_opt', {}) or {}
    g_hf     = cfg.get('gaussian_hf',  {}) or {}
    sl       = cfg.get('slurm', {}) or {}

    job_name = sl.get('job_name') or f"{lig}_ligand_parameters"

    opt_com = Path(g_opt.get('output_com') or f"{lig}_opt.com").name
    hf_com  = Path(g_hf.get('output_com')  or f"{lig}_hf.com").name
    opt_log = Path(opt_com).stem + '.log'
    hf_log  = Path(hf_com).stem  + '.log'

    # -- #SBATCH directives ----------------------------------------------------
    _SKIP = {'job_name', 'modules', 'conda_env', 'cpu'}
    directives = [f'#SBATCH --job-name={job_name}']
    for key, val in sl.items():
        if key in _SKIP or val is None or val == '':
            continue
        directives.append(f'#SBATCH --{key.replace("_", "-")}={val}')

    # -- module load lines -----------------------------------------------------
    modules = sl.get('modules') or []
    module_lines = (f"module load {' '.join(modules)}"
                    if modules else '# (no modules configured)')

    # -- conda activation ------------------------------------------------------
    conda_env = (sl.get('conda_env') or '').strip()
    if conda_env:
        conda_block = (
            f'source "$(conda info --base)/etc/profile.d/conda.sh"\n'
            f'conda activate {conda_env}\n'
        )
    else:
        conda_block = ''

    if config_path is None:
        config_path = Path(proj_root) / 'config.yaml'

    script = _TEMPLATE.format(
        sbatch_directives='\n'.join(directives),
        module_lines=module_lines,
        conda_block=conda_block,
        proj_root=str(proj_root),
        workdir=str(workdir),
        config_path=str(config_path),
        lig=lig,
        opt_com=opt_com,
        opt_log=opt_log,
        hf_com=hf_com,
        hf_log=hf_log,
    )

    output = Path(output)
    output.write_text(script)
    output.chmod(0o755)
    print(f"SLURM script : {output}")
    print(f"  Ligand     : {lig}")
    print(f"  Workdir    : {workdir}")
    for d in directives:
        print(f"  {d.lstrip('#')}")
    if modules:
        print(f"  module load : {', '.join(modules)}")
    if conda_env:
        print(f"  conda env   : {conda_env}")
    print(f"  Workflow    : {opt_com} -> {opt_log} -> {hf_com} -> {hf_log}")
    print(f"\nTo submit:")
    print(f"  sbatch {output}")

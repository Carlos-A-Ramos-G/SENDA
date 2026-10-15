"""Regression tests for the AMBER-path script generation (no AMBER/OpenMM needed)."""

from senda.workflow.run import _param_script

# --- senda-slurm inserts amber_module verbatim ------------------------------

def test_param_script_inserts_amber_module_verbatim():
    slurm = {"amber_module": "module load amber/24", "cpu": {"time": "1:00:00"}}
    script = _param_script(slurm, "config.yaml", "logs")
    assert "\nmodule load amber/24\n" in script
    assert "module load module load" not in script

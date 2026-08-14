"""
senda.qmmm.common.templates

AMBER input file and SLURM script templates for QM/MM string method.
All placeholders use the __KEY__ convention; replaced by fill().
"""


def fill(template: str, **kwargs) -> str:
    """Replace __KEY__ placeholders in template with str(value)."""
    result = template
    for key, value in kwargs.items():
        result = result.replace(f"__{key}__", str(value))
    return result


# ---------------------------------------------------------------------------
# AMBER input: stage 05 (equilibration and restraint-free production)
#
# __NMROPT__  -- either "\n  nmropt   = 1," or "" (no restraints)
# __DISANG__  -- either "&wt type = 'END'/\nDISANG=restr\n/\n" or ""
# __IREST__   -- 0 for fresh start (equil), 1 for restart (restraint-free prod)
# __NTX__     -- 1 for fresh start, 5 for restart with velocities
# ---------------------------------------------------------------------------

STAGE05_IN = """\
QM/MM stage 05
 &cntrl
  temp0    = __TEMP__,
  tempi    = __TEMP__,
  ntt      = 3,
  irest    = __IREST__,
  ntx      = __NTX__,
  ntb      = 1,
  cut      = __QMCUT__,
  gamma_ln = __GAMMA_LN__,
  nstlim   = __NSTLIM__,
  dt       = __DT__,
  ntpr     = __NTPR__,
  ntwx     = __NTWX__,
  ntwr     = __NTWR__,
  ntxo     = 1,
  ifqnt    = 1,__NMROPT__
 /
 &qmmm
  qmmask   = '__QMMASK__',
  qmcharge = __QMCHARGE__,
  qm_theory= '__QM_THEORY__',
  qmcut    = __QMCUT__,
  qm_ewald = 0,
  writepdb = 1,
 /
__DISANG__"""

EQUIL_IN = STAGE05_IN  # kept for any direct imports

# ---------------------------------------------------------------------------
# AMBER input: restrained scan (stage 06, one per node via sed)
# ---------------------------------------------------------------------------

SCAN_IN_TEMPLATE = """\
restrained window __NODE__
 &cntrl
  temp0    = __TEMP__,
  tempi    = __TEMP__,
  ntt      = 3,
  irest    = 0,
  ntx      = 1,
  ntb      = 1,
  cut      = __QMCUT__,
  gamma_ln = __GAMMA_LN__,
  nstlim   = __NSTLIM__,
  dt       = __DT__,
  ntpr     = __NTPR__,
  ntwx     = __NTWX__,
  ntwr     = __NTWR__,
  ntxo     = 1,
  ifqnt    = 1,
  nmropt   = 1,
 /
 &qmmm
  qmmask   = '__QMMASK__',
  qmcharge = __QMCHARGE__,
  qm_theory= '__QM_THEORY__',
  qmcut    = __QMCUT__,
  qm_ewald = 0,
  writepdb = 1,
 /
&wt type = 'END'/
DISANG=restr__NODE__
/
"""

# ---------------------------------------------------------------------------
# AMBER input: string method (stage 07, seed replaced by in.sh)
# ---------------------------------------------------------------------------

STRING_IN = """\
string input file
 &cntrl
  temp0    = __TEMP__,
  tempi    = __TEMP__,
  ntt      = 3,
  irest    = 1,
  ntx      = 5,
  cut      = __QMCUT__,
  gamma_ln = __GAMMA_LN__,
  nstlim   = __NSTLIM__,
  dt       = __DT__,
  ntpr     = __NTPR__,
  ntwx     = __NTWX__,
  ntwr     = __NTWR__,
  ntxo     = 1,
  ifqnt    = 1,
  ig       = @NODE_SEED@,
  asm      = 1,
 /
 &qmmm
  qmmask   = '__QMMASK__',
  qmcharge = __QMCHARGE__,
  qm_theory= '__QM_THEORY__',
  qmcut    = __QMCUT__,
  qm_ewald = 1,
  writepdb = 1,
 /
&wt type = 'END'/
/
 &asm
  preparation_steps = __PREP_STEPS__
  guess_file = 'guess'
  z_bias = .__Z_BIAS__.
  force_constant_d = __FORCE_CONSTANT_D__
"""

# ---------------------------------------------------------------------------
# SLURM: equilibration job (stage 05)
# ---------------------------------------------------------------------------

EQUIL_SLURM = """\
#!/bin/bash
#SBATCH --time=__TIME__
#SBATCH --job-name=__SCHEME___equil
#SBATCH --ntasks=__NTASKS__
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
__ACCOUNT_LINE__
__PARTITION_LINE__

hostname
srun numactl -s
__AMBER_MODULE__

export MPICH_NO_BUFFER_ALIAS_CHECK=1
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

cp __REP_RST7__ 0.rst7
srun --cpu-bind=cores sander.MPI \\
    -O -rem 0 -i in -o 0e.out -c 0.rst7 -r 0e.rst7 \\
    -x 0e.nc -inf 0e.mdinfo -p __PARM__
"""

# ---------------------------------------------------------------------------
# SLURM: restrained scan job (stage 06)
# ---------------------------------------------------------------------------

SCAN_SLURM = """\
#!/bin/bash
#SBATCH --time=__TIME__
#SBATCH --job-name=__SCHEME___scan
#SBATCH --ntasks=__NTASKS__
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
__ACCOUNT_LINE__
__PARTITION_LINE__

hostname
srun numactl -s
__AMBER_MODULE__

export MPICH_NO_BUFFER_ALIAS_CHECK=1
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

NODES=__N_NODES__

if [ -f ../05_QMMM_restraint_free/prod.rst7 ]; then
  cp ../05_QMMM_restraint_free/prod.rst7 0.rst7
else
  cp ../05_QMMM_equilibration/0e.rst7 0.rst7
fi

for i in $(seq 1 $NODES); do
  sed -e "s/__NODE__/$i/g" in_template > in
  cat restr0 >> restr${i}
  srun --cpu-bind=cores sander.MPI \\
      -O -rem 0 -i in -o ${i}.out -c $((i-1)).rst7 -r ${i}.rst7 \\
      -x ${i}.nc -inf ${i}.mdinfo -p __PARM__
  if [ $? -ne 0 ]; then
    echo "ERROR: sander.MPI failed at node $i -- aborting." >&2
    exit 1
  fi
done

bash center.sh
"""

# ---------------------------------------------------------------------------
# cpptraj: centering script run at end of scan job
# ---------------------------------------------------------------------------

CENTER_SH = """\
#!/bin/bash
# Center scan restart files on the protein before string method.

PARM=__PARM__
N=__N_NODES__
MASK=":1-__PROTEIN_LAST_RES__"

for i in $(seq 1 $N); do
  INPUT="${i}.rst7"
  OUTPUT="${i}_centred.rst7"

  if [ ! -f "$INPUT" ]; then
    echo "WARNING: $INPUT not found, skipping."
    continue
  fi

  cpptraj <<EOF
parm $PARM
trajin $INPUT
autoimage
center $MASK mass origin
image origin center
trajout $OUTPUT restart vel
go
quit
EOF

  if [ $? -ne 0 ]; then
    echo "ERROR: cpptraj failed on $INPUT -- aborting." >&2
    exit 1
  fi
done

echo "Centering complete: ${N} files written as *_centred.rst7"
"""

# ---------------------------------------------------------------------------
# Bash: string node setup script (stage 07)
# ---------------------------------------------------------------------------

STRING_IN_SH = """\
#!/bin/bash
# Generate per-node input files and groupfile for sander.MPI -ng.

if [ -f string.groupfile ]; then
  rm string.groupfile
fi

NODES=__N_NODES__
PARM=__PARM_H10__
REACT=../06_QMMM_scan/1_centred.rst7
PROD=../06_QMMM_scan/${NODES}_centred.rst7
SEED=__SEED__

for i in $(seq 1 $NODES); do
  # generate sander input for node $i
  sed "s/@NODE_SEED@/$((SEED + i))/g" in > ${i}.in

  # use the reactant structure for the first half of the nodes, product for the rest
  if [ $i -le $((NODES / 2)) ]; then
    crd=$REACT
  else
    crd=$PROD
  fi

  # write out the sander arguments for node $i to the groupfile
  echo "-O -rem 0 -i ${i}.in -o ${i}.out -c ${crd} -r ${i}.rst7 -x ${i}.nc -inf ${i}.mdinfo -p $PARM" >> string.groupfile
done
"""

# ---------------------------------------------------------------------------
# SLURM: string method job (stage 07)
# ---------------------------------------------------------------------------

STRING_SLURM = """\
#!/bin/bash
#SBATCH --time=__TIME__
#SBATCH --job-name=__SCHEME___string
#SBATCH --ntasks=__NTASKS_STRING__
__ACCOUNT_LINE__
__PARTITION_LINE__

hostname
srun numactl -s
__AMBER_MODULE__

export MPICH_NO_BUFFER_ALIAS_CHECK=1

mkdir -p results
bash in.sh __N_NODES__
srun --cpu-bind=cores sander.MPI -ng __N_NODES__ -groupfile string.groupfile
"""

# ---------------------------------------------------------------------------
# SLURM: restraint-free production job (stage 05_QMMM_restraint_free)
# ---------------------------------------------------------------------------

PROD_SLURM = """\
#!/bin/bash
#SBATCH --time=__TIME__
#SBATCH --job-name=__SCHEME___prod
#SBATCH --ntasks=__NTASKS__
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
__ACCOUNT_LINE__
__PARTITION_LINE__

hostname
srun numactl -s
__AMBER_MODULE__

export MPICH_NO_BUFFER_ALIAS_CHECK=1
export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

cp ../05_QMMM_equilibration/0e.rst7 0.rst7
srun --cpu-bind=cores sander.MPI \\
    -O -rem 0 -i in -o prod.out -c 0.rst7 -r prod.rst7 \\
    -x prod.nc -inf prod.mdinfo -p __PARM__
"""

#!/bin/bash
# Build the OpenMM replica from the AMBER run's 00_prep (no tleap: same parm7,
# coordinates and restraints) with the protocol of the AMBER run's config.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need python; need cpptraj
python -c "import openmm, senda" || { echo "ERROR: python needs openmm and senda"; exit 1; }

restr_arg=()
if [ -f "${RESTR}" ]; then
    restr_arg=(--restr "${RESTR}")
else
    log "no restr file at ${RESTR} -- assuming the AMBER run had no NMR restraints"
fi

$COMPARE prepare-openmm \
    --prep "${PREP}" ${restr_arg[@]+"${restr_arg[@]}"} --config "${CONFIG}" \
    --out "${OPENMM_RUN}" --prod-ns "${PROD_NS}" --chunk-ns "${CHUNK_NS}" \
    --platform "${PLATFORM}" --seed "${SEED}" --openmm-command "${OMM_COMMAND}" "$@"

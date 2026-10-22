#!/bin/bash
# Compare the OpenMM run with the AMBER run over the same production window
# (first PROD_NS ns of each). Writes report/summary.md and report/comparison.png.
source "$(dirname "${BASH_SOURCE[0]}")/_common.sh"
need python; need cpptraj

$COMPARE compare --amber "${AMBER_REPLICA}" --openmm "${OPENMM_RUN}" \
    --out "${REPORT}" --prod-ns "${PROD_NS}" "$@"

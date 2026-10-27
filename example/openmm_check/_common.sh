# Sourced by the numbered scripts: loads settings.sh, the environment and paths.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"
source ./settings.sh
eval "${ENV_SETUP}"      # before set -u: conda's activate scripts use unset variables
set -u

SENDA_ROOT="$(cd "$HERE/../.." && pwd)"
COMPARE="python ${SENDA_ROOT}/scripts/compare_engines.py"
[ -n "${RESTR}" ] || RESTR="${AMBER_REPLICA}/restr"

need() { command -v "$1" >/dev/null || { echo "ERROR: '$1' not on PATH (see ENV_SETUP in settings.sh)"; exit 1; }; }
log()  { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

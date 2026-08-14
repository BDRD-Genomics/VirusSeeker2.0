#!/usr/bin/env bash
set -euo pipefail

REPEATMASKER_ENV="/opt/conda/envs/repeatmasker"

# FamDB uses Python dependencies installed in the RepeatMasker environment.
# Put this environment first even when VirusSeeker activated the vs environment.
export PATH="${REPEATMASKER_ENV}/bin:${PATH}"
export PYTHONNOUSERSITE=1

exec "${REPEATMASKER_ENV}/bin/RepeatMasker" "$@"

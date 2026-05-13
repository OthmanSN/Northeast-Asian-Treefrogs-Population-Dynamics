#!/usr/bin/env bash
set -euo pipefail

# Generate ngsPSMC input from one BAM file using ANGSD.
# Usage:
#   bash scripts/04_generate_ngspsmc_input.sh /path/to/angsd /path/to/sample.bam results/ngspsmc/sample

ANGSD=${1:?"Path to ANGSD executable is required"}
BAM=${2:?"Input BAM file is required"}
OUT_PREFIX=${3:?"Output prefix is required"}

mkdir -p "$(dirname "$OUT_PREFIX")"

"$ANGSD" \
  -i "$BAM" \
  -dopsmc 1 \
  -out "$OUT_PREFIX" \
  -gl 1 \
  -minq 20 \
  -minmapq 30

echo "Done. ngsPSMC input prefix: $OUT_PREFIX"

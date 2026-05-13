#!/usr/bin/env bash
set -euo pipefail

# Run ABBA-BABA tests with Dsuite Dtrios, optional Dinvestigate, and f-branch statistics.
# Usage:
#   bash scripts/02_run_dsuite_abba_baba.sh config/config.sh

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$OUTDIR"/{dsuite,logs}

SORTED_VCF="$OUTDIR/dsuite/populations_snps_final.sorted.vcf.gz"
PREFIX="$OUTDIR/dsuite/dstats"
FBRANCH_OUT="$OUTDIR/dsuite/fbranch_out.txt"

# 1. Sort and index VCF for Dsuite.
"$BCFTOOLS" sort "$RAW_VCF" -Oz -o "$SORTED_VCF"
"$TABIX" -f -p vcf "$SORTED_VCF"

# 2. Run Dtrios with a species tree.
"$DSUITE" Dtrios \
  -t "$SPECIES_TREE" \
  -n dstats \
  "$SORTED_VCF" \
  "$SPECIES_SETS" \
  > "$OUTDIR/logs/dsuite_dtrios.log" 2>&1

# Dsuite writes output files to the working directory depending on version and prefix.
# Move common outputs into the dsuite result folder when present.
for file in dstats_BBAA.txt dstats_tree.txt dstats_combine.txt; do
  [[ -f "$file" ]] && mv "$file" "$OUTDIR/dsuite/"
done

# 3. Compute f-branch statistics if Dtrios tree output exists.
if [[ -f "$OUTDIR/dsuite/dstats_tree.txt" ]]; then
  "$DSUITE" Fbranch "$SPECIES_TREE" "$OUTDIR/dsuite/dstats_tree.txt" > "$FBRANCH_OUT"
else
  echo "Warning: dstats_tree.txt not found; skipping Fbranch." >&2
fi

# 4. Optional local-window introgression scan using Dinvestigate.
if [[ -f "$TEST_TRIOS" ]]; then
  "$DSUITE" Dinvestigate \
    "$SORTED_VCF" \
    "$SPECIES_SETS" \
    "$TEST_TRIOS" \
    -w 100000,50000 \
    > "$OUTDIR/logs/dsuite_dinvestigate.log" 2>&1
fi

echo "Done. Dsuite outputs are in: $OUTDIR/dsuite"

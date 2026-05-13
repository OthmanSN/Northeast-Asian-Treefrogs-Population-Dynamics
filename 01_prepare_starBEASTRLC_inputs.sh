#!/usr/bin/env bash
set -euo pipefail

# Prepare RADseq SNP datasets for SNP-based divergence dating and StarBEASTRLC input.
# Usage:
#   bash scripts/01_prepare_starBEESTRLC_inputs.sh config/config.sh

CONFIG=${1:-config/config.sh}
source "$CONFIG"

mkdir -p "$OUTDIR"/{filtered_vcfs,loci_vcfs,loci_nexus,logs}

FILTERED_VCF="$OUTDIR/filtered_vcfs/snps_biallelic_missing30.vcf.gz"
FILTERED_PLAIN_VCF="$OUTDIR/filtered_vcfs/snps_biallelic_missing30.vcf"
WINDOWS="$OUTDIR/loci_vcfs/windows_min15.txt"

# 1. Keep biallelic SNPs and remove SNPs with >30% missing data.
"$BCFTOOLS" view \
  -i 'F_MISSING<0.3 && N_ALT=1' \
  -v snps \
  "$RAW_VCF" \
  -Oz -o "$FILTERED_VCF"
"$BCFTOOLS" index -f "$FILTERED_VCF"

# 2. Basic QC summaries.
{
  echo "Filtered VCF: $FILTERED_VCF"
  echo -n "Number of SNPs: "
  "$BCFTOOLS" view -v snps -H "$FILTERED_VCF" | wc -l
  echo -n "Number of indels: "
  "$BCFTOOLS" view -v indels -H "$FILTERED_VCF" | wc -l
} > "$OUTDIR/logs/filtered_vcf_summary.txt"

# 3. Convert compressed VCF to plain VCF for tools/scripts requiring uncompressed input.
"$BCFTOOLS" view -Ov -o "$FILTERED_PLAIN_VCF" "$FILTERED_VCF"

# 4. Identify RAD loci/windows with at least 15 variants in 500 bp windows.
"$VCFTOOLS" --vcf "$FILTERED_PLAIN_VCF" \
  --window-pi 500 \
  --out "$OUTDIR/loci_vcfs/windowed_output"

awk 'NR>1 && $4>=15 {print $1, $2, $3}' \
  "$OUTDIR/loci_vcfs/windowed_output.windowed.pi" > "$WINDOWS"

# 5. Split VCF by selected windows/loci.
while read -r chrom start end; do
  outfile="$OUTDIR/loci_vcfs/locus_${chrom}_${start}_${end}"
  "$VCFTOOLS" --vcf "$FILTERED_PLAIN_VCF" \
    --chr "$chrom" \
    --from-bp "$start" \
    --to-bp "$end" \
    --recode --out "$outfile"
done < "$WINDOWS"

# 6. Convert each locus VCF to NEXUS for downstream dating workflows.
for vcf in "$OUTDIR"/loci_vcfs/*.recode.vcf; do
  "$PYTHON" "$VCF2PHYLIP" \
    -i "$vcf" \
    -n \
    --output-folder "$OUTDIR/loci_nexus"
done

echo "Done. StarBEASTRLC-ready locus files are in: $OUTDIR/loci_nexus"

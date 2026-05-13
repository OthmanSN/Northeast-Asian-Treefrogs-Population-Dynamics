#!/usr/bin/env bash
# Copy this file to config.sh and edit paths for your local machine.

# Project directories
PROJECT_DIR="/path/to/project"
RAW_VCF="${PROJECT_DIR}/data/populations_snps_final.vcf.gz"
OUTDIR="${PROJECT_DIR}/results"

# External tools
BCFTOOLS="bcftools"
VCFTOOLS="vcftools"
TABIX="tabix"
BGZIP="bgzip"
DSUITE="Dsuite"
PYTHON="python3"

# Helper scripts
VCF2PHYLIP="/path/to/vcf2phylip.py"
DSUITE_UTILS_DIR="/path/to/Dsuite/utils"

# Analysis inputs
SPECIES_SETS="${PROJECT_DIR}/metadata/species_sets.txt"
SPECIES_TREE="${PROJECT_DIR}/metadata/species_tree.nwk"
TEST_TRIOS="${PROJECT_DIR}/metadata/test_trios.txt"

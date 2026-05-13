# Dryophytes RADseq downstream workflows

This repository contains cleaned command-line workflows used for RADseq downstream analyses, including SNP filtering for divergence dating, preparation of StarBEASTRLC locus inputs, ABBA–BABA introgression tests with Dsuite, f-branch plotting, and ngsPSMC input generation.

## Repository structure

```text
.
├── config/
│   └── config.template.sh
├── scripts/
│   ├── 01_prepare_starBEESTRLC_inputs.sh
│   ├── 02_run_dsuite_abba_baba.sh
│   ├── 03_plot_fbranch_heatmap.py
│   └── 04_generate_ngspsmc_input.sh
└── README.md
```

## Requirements

Install or make available in your PATH:

- `bcftools`
- `vcftools`
- `tabix` and `bgzip` from HTSlib
- `Dsuite`
- `python >= 3.8`
- Python packages: `pandas`, `matplotlib`, `seaborn`
- `vcf2phylip.py`
- Optional: `ANGSD` for ngsPSMC input generation

Example Python setup:

```bash
pip install pandas matplotlib seaborn
```

## Configuration

Copy the template and edit the paths for your local machine:

```bash
cp config/config.template.sh config/config.sh
nano config/config.sh
```

The file `config/config.sh` should contain local paths to the input VCF, output directory, Dsuite, the species tree, species-set file, and helper scripts. Do not commit `config/config.sh` if it contains local or private paths.

## 1. SNP filtering and StarBEASTRLC locus preparation

```bash
bash scripts/01_prepare_starBEESTRLC_inputs.sh config/config.sh
```

This script performs the following steps:

1. Filters to biallelic SNPs with less than 30% missing data.
2. Indexes the filtered VCF and writes a short QC summary.
3. Converts the VCF to uncompressed format when required by downstream tools.
4. Identifies 500 bp windows with at least 15 variants.
5. Splits the VCF into locus-level files.
6. Converts locus-level VCFs to NEXUS format for dating workflows.

## 2. ABBA–BABA and f-branch analyses with Dsuite

```bash
bash scripts/02_run_dsuite_abba_baba.sh config/config.sh
```

This script sorts and indexes the VCF, runs `Dsuite Dtrios`, computes f-branch statistics using `Dsuite Fbranch` when the tree output is available, and optionally runs `Dsuite Dinvestigate` if a trio file is provided.

## 3. Plot f-branch heatmap

```bash
python scripts/03_plot_fbranch_heatmap.py \
  --input results/dsuite/fbranch_out.txt \
  --output results/dsuite/fbranch_heatmap.png
```

## 4. Generate ngsPSMC input with ANGSD

```bash
bash scripts/04_generate_ngspsmc_input.sh \
  /path/to/angsd \
  /path/to/sample.bam \
  results/ngspsmc/sample_name
```

## Notes

- These scripts are templates for reproducible workflows. Before running, check all paths, file names, taxon names, and Dsuite input formats.
- The original workflow used several exploratory thresholds. This cleaned version keeps the main filtering threshold at 30% missing data and a minimum of 15 SNPs per 500 bp locus window. Adjust these values if the manuscript requires different filtering criteria.
- Keep raw data, large VCF/BAM files, and private local paths out of GitHub. Use `.gitignore` for large outputs and machine-specific configuration files.

## Suggested citation

If these scripts support a publication, cite the relevant software used in the workflow, including BCFtools, VCFtools, Dsuite, ANGSD, and vcf2phylip, alongside the manuscript or data repository citation.

# SNP and Sample Quality Control (QC)

## Purpose
This tool performs quality control filtering on GWAS genotype data.  
It removes SNPs and samples that do not meet thresholds for:
- Minor Allele Frequency (MAF)
- SNP missingness
- Sample missingness

Ensures only high-quality variants and samples are used for downstream analysis.

## Inputs
- **Genotypes** (tabular): Matrix of genotypes (samples × SNPs).
- **SNP Annotation** (tabular): Information about each SNP (ID, chromosome, position, MAF, etc.).
- **Phenotypes** (tabular, optional): Trait and covariate data for each sample.
- **MAF Threshold** (float): Minimum allowed minor allele frequency.
- **SNP Missingness Threshold** (float): Maximum allowed missing genotype rate for SNPs.
- **Sample Missingness Threshold** (float): Maximum allowed missing genotype rate for samples.
- **Missing Value Code** (string): Code representing missing values in the genotype matrix (e.g., `NA`).

## Outputs
- **Filtered Genotypes** (tabular): QC-passed genotype matrix.
- **Filtered SNP Annotation** (tabular): QC-passed SNP annotation file.
- **Filtered Samples** (tabular): QC-passed sample list (or phenotype file if provided).
- **QC Report** (txt): Summary of filtering steps and statistics before/after QC.

## Test
```bash
python tools/snp_sample_qc/snp_sample_qc.py \
  --genotypes gwas_data/genotypes.tsv \
  --snp_annotation gwas_data/snp_annotation.tsv \
  --phenotypes gwas_data/phenotypes_covariates.tsv \
  --maf_thresh 0.2 \
  --snp_miss_thresh 0.05 \
  --out_genotypes outputs/snp_sample_qc/geno_filtered.tsv \
  --out_snp_annotation outputs/snp_sample_qc/snp_annot_filtered.tsv \
  --out_samples outputs/snp_sample_qc/samples_filtered.tsv \
  --qc_report outputs/snp_sample_qc/qc_report.txt

# GWAS Association Testing

## Purpose
Performs genome-wide association testing for each SNP in the filtered genotype data against the provided phenotype.  

## Inputs
- **Genotypes** (tabular): QC-passed genotype matrix.
- **Phenotypes** (tabular): File with sample IDs and phenotype values.
- **SNP Annotation** (tabular): QC-passed SNP metadata.
- **Covariates** (comma-separated list, optional): Covariates to adjust for.
- **Minimum MAF** (float): Minimum allowed minor allele frequency.
- **Missing Value Code** (string): Code representing missing values in the genotype matrix (e.g., `NA`).

## Outputs
- **GWAS Results** (tabular): One row per SNP with:
  - SNP ID
  - Chromosome
  - Position
  - Beta
  - Standard error
  - p-value
  - MAF

## Test
```bash
python tools/association_test/association_test.py \
  --genotypes outputs/snp_sample_qc/geno_filtered.tsv \
  --snp_annotation outputs/snp_sample_qc/snp_annot_filtered.tsv \
  --phenotypes gwas_data/phenotypes_covariates.tsv \
  --covariates age,sex \
  --out outputs/association_test/assoc_results.tsv


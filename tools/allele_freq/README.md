# Allele Frequency Calculation

## Purpose
- Calculates allele frequencies for each SNP in the genotype dataset.
- Outputs per-SNP counts of genotypes (0,1,2), n_non_missing, MAF.
- If n_non_missing < 10, MAF will be set to NA.
- Helps to identify rare and common variants for filtering and interpretation.

## Inputs
- **Genotypes** (tabular): QC-passed genotype matrix.
- **SNP Annotation** (tabular): SNP metadata including ID and alleles.
- **Missing Value Code** (string): Code representing missing values in the genotype matrix (e.g., `NA`).

## Outputs
- **Allele Frequency Table** (tabular): For each SNP, includes:
  - SNP ID
  - Chromosome
  - Position
  - Allele counts
  - Allele frequencies (MAF, major/minor alleles)

## Test
```bash
python tools/allele_freq/allele_freq.py \
  --genotypes outputs/snp_sample_qc/geno_filtered.tsv \
  --snp_annotation outputs/snp_sample_qc/snp_annot_filtered.tsv \
  --out outputs/allele_freq/allele_freqs.tsv

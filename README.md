# GWAS Galaxy Tools

A collection of **five Galaxy Python tools** for a **Genome-Wide Association Study (GWAS)** analysis.  
Designed for use with genotype/phenotype datasets.

## Tools Overview
1. **SNP & Sample QC** — Filter genotype data by MAF and missingness thresholds.
2. **Allele Frequency Calculator** — Compute per-SNP genotype counts and minor allele frequency (MAF).
3. **Association Test** — Run per-SNP logistic regression between genotype and binary phenotype, with optional covariates.
4. **Manhattan Plot Generator** — Create a Manhattan plot from association results and output top hits.
5. **Windowed LD Calculator** — Compute LD (r²) for SNPs within a genomic window around a focal SNP.

---

## Project Structure

GWAS-Galaxy-Tools/
```plaintext
├── tools/
│ ├── snp_sample_qc/ 
│ ├── allele_freq/
│ ├── association_test/ 
│ ├── manhattan/ 
│ └── window_ld/ 
│
├── gwas_data/
│ ├── genotypes.tsv
│ ├── phenotypes_covariates.tsv
│ ├── snp_annotation.tsv
│
├── outputs/ 
│ ├── snp_sample_qc/ 
│ ├── allele_freq/
│ ├── association_test/ 
│ ├── manhattan/ 
│ └── window_ld/
├── generate_gwas_data.py 
├── run_log.txt
└── README.md
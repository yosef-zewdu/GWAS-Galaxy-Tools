# Local Linkage Disequilibrium (LD) Calculation

## Purpose
Calculates pairwise LD for SNPs within a genomic window centered on a focal SNP.  
Generates both a matrix and heatmap.

## Inputs
- **Genotypes** (tabular): QC-passed genotype matrix.
- **SNP Annotation** (tabular): QC-passed SNP annotation file.
- **Focal SNP** (string): SNP ID to center the analysis.
- **Window (kb)** (integer): Size of genomic window (± kb).

## Outputs
- **LD Matrix** (tabular): Matrix of LD values (r²) between SNPs.
- **LD Heatmap** (png): Heatmap of LD matrix.

## Test
```bash
python tools/window_ld/window_ld.py \
  --genotypes outputs/snp_sample_qc/geno_filtered.tsv \
  --snp_annotation outputs/snp_sample_qc/snp_annot_filtered.tsv \
  --focal_snp rs100696 \
  --window_kb 250 \
  --out_matrix outputs/window_ld/ld_matrix.tsv \
  --out_png outputs/window_ld/ld_heatmap.png


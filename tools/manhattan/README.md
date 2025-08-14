# Manhattan Plot Visualization

## Purpose
Creates a Manhattan plot from GWAS association results, showing -log10(p-values) across the genome.

## Inputs
- **Association Results** (tabular): GWAS results table with columns: `SNP_ID`, `CHROM`, `pos`, `P`, `MAF` etc.
- **SNP Annotation** (tabular): QC-passed SNP metadata.
- **P-value Threshold** (float): Significance threshold (default 5e-8).
- **Top N hits** (integer):  Number of top SNPs to report in the output table (default 20).

## Outputs
- **Manhattan Plot** (png): Visualization of SNP associations.
- **Top hits** (tabular): Table of the top SNPs hit by p-value. 

## Test
```bash
python tools/manhattan/manhattan.py \
  --assoc outputs/association_test/assoc_results.tsv \
  --snp_annotation outputs/snp_sample_qc/snp_annot_filtered.tsv \
  --p_threshold 1e-4 \
  --out_png outputs/manhattan/manhattan.png \
  --out_top outputs/manhattan/top_hits.tsv
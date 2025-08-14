import argparse
import pandas as pd
import numpy as np
import sys

def compute_maf_from_genotypes(col):
    # col contains 0/1/2 or NaN; allele frequency of alt = sum(genotypes)/(2*n_nonmissing)
    nonmiss = col.dropna()
    n = nonmiss.shape[0]
    if n == 0:
        return np.nan
    alt_af = nonmiss.sum() / (2.0 * n)
    maf = min(alt_af, 1.0 - alt_af)
    return maf

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genotypes", required=True, help="Input text file")
    p.add_argument("--snp_annotation", required=True)
    p.add_argument("--phenotypes", default=None)
    p.add_argument("--maf_thresh", type=float, default=0.01)
    p.add_argument("--snp_miss_thresh", type=float, default=0.05)
    p.add_argument("--sample_miss_thresh", type=float, default=0.1)
    p.add_argument("--missing_code", default="NA")
    p.add_argument("--out_genotypes", required=True)
    p.add_argument("--out_snp_annotation", required=True)
    p.add_argument("--out_samples", required=True)
    p.add_argument("--qc_report", required=True)
    args = p.parse_args()

    # Read genotype matrix (samples x SNPs)
    geno = pd.read_csv(args.genotypes, sep="\t", index_col=0)
    # Replace specified missing code with NaN
    geno = geno.replace(args.missing_code, np.nan).astype(float)

    snp_ann = pd.read_csv(args.snp_annotation, sep="\t")
    phenos = pd.read_csv(args.phenotypes, sep="\t") if args.phenotypes else None

    # SNP missingness: fraction missing per SNP (column)
    snp_missingness = geno.isna().mean(axis=0)
    # Sample missingness: fraction missing per sample (row)
    sample_missingness = geno.isna().mean(axis=1)

    # Compute MAF per SNP
    maf_series = geno.apply(compute_maf_from_genotypes, axis=0)

    # Filter SNPs by missingness and MAF
    snp_keep_mask = (snp_missingness <= args.snp_miss_thresh) & (maf_series >= args.maf_thresh)
    kept_snps = snp_keep_mask[snp_keep_mask].index.tolist()

    # Filter samples by missingness
    sample_keep_mask = (sample_missingness <= args.sample_miss_thresh)
    kept_samples = sample_keep_mask[sample_keep_mask].index.tolist()

    # Apply filters
    geno_filtered = geno.loc[kept_samples, kept_snps]
    snp_ann_filtered = snp_ann[snp_ann['snp_id'].isin(kept_snps)].copy()
    phenos_filtered = None
    if phenos is not None:
        # Expect sample_id column
        if 'sample_id' not in phenos.columns:
            sys.exit("phenotypes file must have 'sample_id' column")
        phenos_filtered = phenos[phenos['sample_id'].isin(kept_samples)]

    # QC report
    report_lines = []
    report_lines.append(f"Initial SNP count: {geno.shape[1]}")
    report_lines.append(f"Initial sample count: {geno.shape[0]}")
    report_lines.append(f"SNP missingness threshold: {args.snp_miss_thresh}")
    report_lines.append(f"Sample missingness threshold: {args.sample_miss_thresh}")
    report_lines.append(f"MAF threshold: {args.maf_thresh}")
    report_lines.append("")
    report_lines.append(f"SNPs passing filters: {len(kept_snps)}")
    report_lines.append(f"Samples passing filters: {len(kept_samples)}")
    report_lines.append("")
    report_lines.append("Top 10 SNPs removed by missingness (snp -> missing_frac):")
    report_lines += [f"{s} -> {snp_missingness[s]:.4f}" for s in snp_missingness.sort_values(ascending=False).head(10).index]
    report_lines.append("")
    report_lines.append("Top 10 SNPs removed by low MAF (snp -> maf):")
    report_lines += [f"{s} -> {maf_series[s]}" for s in maf_series.sort_values().head(10).index]

    # Save outputs
    geno_filtered.to_csv(args.out_genotypes, sep="\t", index=True)
    snp_ann_filtered.to_csv(args.out_snp_annotation, sep="\t", index=False)
    if phenos_filtered is not None:
        phenos_filtered.to_csv(args.out_samples, sep="\t", index=False)
    else:
        # Output sample list (one per line)
        pd.DataFrame({'sample_id': kept_samples}).to_csv(args.out_samples, sep="\t", index=False)
    with open(args.qc_report, "w") as fh:
        fh.write("\n".join(report_lines))

    print("QC done. Outputs written.")

if __name__ == "__main__":
    main()

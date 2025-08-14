import argparse
import pandas as pd
import numpy as np

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genotypes", required=True)
    p.add_argument("--snp_annotation", default=None)
    p.add_argument("--sample_subset", default=None)
    p.add_argument("--missing_code", default="NA")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    geno = pd.read_csv(args.genotypes, sep="\t", index_col=0)
    geno = geno.replace(args.missing_code, np.nan).astype(float)

    if args.sample_subset:
        subset = pd.read_csv(args.sample_subset, sep="\t", header=None)
        samples = subset.iloc[:,0].astype(str).tolist()
        geno = geno.loc[geno.index.intersection(samples), :]

    records = []
    for snp in geno.columns:
        col = geno[snp].dropna()
        n_non = col.shape[0]
        c0 = (col == 0).sum()
        c1 = (col == 1).sum()
        c2 = (col == 2).sum()
        if n_non < 10:
            maf = np.nan
        else:
            alt_af = col.sum() / (2.0 * n_non)
            maf = min(alt_af, 1.0 - alt_af)
        records.append((snp, c0, c1, c2, n_non, maf))

    df = pd.DataFrame(records, columns=['snp_id', 'count_0', 'count_1', 'count_2', 'n_non_missing', 'MAF'])
    if args.snp_annotation:
        ann = pd.read_csv(args.snp_annotation, sep="\t")
        df = df.merge(ann, left_on='snp_id', right_on='snp_id', how='left')
        # Reorder columns
        cols = ['snp_id', 'chrom', 'pos', 'count_0','count_1','count_2','n_non_missing','MAF']
        df = df[cols]
    df.to_csv(args.out, sep="\t", index=False)
    print("Allele frequency table written to", args.out)

if __name__ == "__main__":
    main()

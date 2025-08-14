import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

def pairwise_r2(x, y):
    # x and y are pandas Series with possible NaNs; compute r^2 using pairwise complete observations
    df = pd.concat([x, y], axis=1).dropna()
    if df.shape[0] < 4:
        return np.nan
    a = df.iloc[:,0].values
    b = df.iloc[:,1].values
    # if constant, correlation undefined
    if np.std(a) == 0 or np.std(b) == 0:
        return np.nan
    r = np.corrcoef(a, b)[0,1]
    return r*r

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genotypes", required=True)
    p.add_argument("--snp_annotation", required=True)
    p.add_argument("--focal_snp", default=None)
    p.add_argument("--chrom", default=None)
    p.add_argument("--center_pos", type=int, default=None)
    p.add_argument("--window_kb", type=int, default=250)
    p.add_argument("--min_maf", type=float, default=0.01)
    p.add_argument("--missing_code", default="NA")
    p.add_argument("--out_matrix", required=True)
    p.add_argument("--out_png", required=True)
    args = p.parse_args()

    geno = pd.read_csv(args.genotypes, sep="\t", index_col=0)
    geno = geno.replace(args.missing_code, np.nan).astype(float)
    ann = pd.read_csv(args.snp_annotation, sep="\t")

    if args.focal_snp:
        if args.focal_snp not in ann['snp_id'].values:
            raise ValueError("focal_snp not found in snp_annotation")
        focal = ann[ann['snp_id'] == args.focal_snp].iloc[0]
        chrom = focal['chrom']
        center_pos = int(focal['pos'])
    else:
        if args.chrom is None or args.center_pos is None:
            raise ValueError("Either focal_snp or chrom+center_pos must be provided")
        chrom = args.chrom
        center_pos = args.center_pos

    window_bp = args.window_kb * 1000
    ann_sel = ann[(ann['chrom'].astype(str) == str(chrom)) & (ann['pos'] >= center_pos - window_bp) & (ann['pos'] <= center_pos + window_bp)].copy()
    if ann_sel.shape[0] == 0:
        raise ValueError("No SNPs in window")
    # Filter by MAF
    # compute MAF from genotype matrix if possible, else use annotation column
    maf_from_ann = 'maf' in ann_sel.columns
    if maf_from_ann:
        ann_sel = ann_sel[ann_sel['maf'] >= args.min_maf]
    else:
        # compute from genotypes
        mafs = {}
        for s in ann_sel['snp_id']:
            if s not in geno.columns:
                mafs[s] = np.nan
            else:
                col = geno[s].dropna()
                n = col.shape[0]
                if n < 4:
                    mafs[s] = np.nan
                else:
                    alt_af = col.sum() / (2.0 * n)
                    mafs[s] = min(alt_af, 1-alt_af)
        ann_sel['maf'] = ann_sel['snp_id'].map(mafs)
        ann_sel = ann_sel[ann_sel['maf'] >= args.min_maf]

    sel_snps = ann_sel.sort_values('pos')['snp_id'].tolist()
    # build LD matrix
    m = len(sel_snps)
    ld = pd.DataFrame(np.full((m,m), np.nan), index=sel_snps, columns=sel_snps)
    for i in range(m):
        for j in range(i, m):
            s1 = sel_snps[i]; s2 = sel_snps[j]
            if s1 not in geno.columns or s2 not in geno.columns:
                val = np.nan
            else:
                val = pairwise_r2(geno[s1], geno[s2])
            ld.iloc[i,j] = val
            ld.iloc[j,i] = val

    ld.to_csv(args.out_matrix, sep="\t", index=True)

    # Heatmap
    fig, ax = plt.subplots(figsize=(6,6))
    im = ax.imshow(ld.values.astype(float), interpolation='nearest', vmin=0, vmax=1)
    ax.set_xticks(range(m))
    ax.set_yticks(range(m))
    ax.set_xticklabels([s for s in sel_snps], rotation=90, fontsize=6)
    ax.set_yticklabels([s for s in sel_snps], fontsize=6)
    fig.colorbar(im, ax=ax, label='r^2')
    plt.tight_layout()
    plt.savefig(args.out_png,format='png', dpi=150)
    plt.close(fig)
    print("LD matrix written to", args.out_matrix)
    print("LD heatmap written to", args.out_png)

if __name__ == "__main__":
    main()

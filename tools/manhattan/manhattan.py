import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--assoc", required=True)
    p.add_argument("--snp_annotation", default=None)
    p.add_argument("--p_threshold", type=float, default=5e-8)
    p.add_argument("--top_n", type=int, default=20)
    p.add_argument("--out_png", required=True)
    p.add_argument("--out_top", required=True)
    args = p.parse_args()

    df = pd.read_csv(args.assoc, sep="\t")
    if 'p' not in df.columns:
        raise ValueError("association file must contain 'p' column")
    if args.snp_annotation and ('chrom' not in df.columns or 'pos' not in df.columns):
        # try merging annotation
        ann = pd.read_csv(args.snp_annotation, sep="\t")
        df = df.merge(ann, on='snp_id', how='left')

    # Prepare plot ordering: sort by chrom and pos
    df['chrom'] = df['chrom'].astype(str)
    # Make chromosome order consistent (1..22 then others)
    try:
        chrom_order = sorted(df['chrom'].unique(), key=lambda x: int(x) if x.isdigit() else 1e6 + hash(x))
    except:
        chrom_order = sorted(df['chrom'].unique())
    df['chrom'] = pd.Categorical(df['chrom'], categories=chrom_order, ordered=True)
    df = df.sort_values(['chrom','pos'])

    # Create a cumulative position for plotting
    chrom_sizes = df.groupby('chrom', observed=False)['pos'].max().to_dict()
    cumulative = {}
    cum = 0
    for c in chrom_order:
        cumulative[c] = cum
        if c in chrom_sizes:
            cum += chrom_sizes[c] + 1e6  # small gap between chromosomes

    df['cpos'] = df.apply(lambda r: r['pos'] + cumulative.get(r['chrom'],0), axis=1)
    df['logp'] = -np.log10(df['p'].replace(0, 1e-300))

    # Plot
    fig, ax = plt.subplots(figsize=(12,5))
    colors = ['#1f77b4', '#ff7f0e']
    for i, (chrom, group) in enumerate(df.groupby('chrom', observed=False)):
        ax.scatter(group['cpos'], group['logp'], s=6, color=colors[i % 2], label=str(chrom))

    ax.axhline(-np.log10(args.p_threshold), color='red', linestyle='--', linewidth=1)
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('-log10(p)')
    # tick positions at center of each chromosome
    ticks = []
    labels = []
    for chrom in chrom_order:
        g = df[df['chrom']==chrom]
        if g.shape[0]==0: continue
        ticks.append(g['cpos'].median())
        labels.append(str(chrom))
    ax.set_xticks(ticks)
    ax.set_xticklabels(labels, rotation=45)
    ax.set_title('Manhattan plot')
    plt.tight_layout()
    plt.savefig(args.out_png, format='png', dpi=150)
    plt.close(fig)

    # Top hits
    top = df.nsmallest(args.top_n, 'p')[['snp_id','chrom','pos','p','logp','MAF' if 'MAF' in df.columns else None]]
    # drop None column if present
    top = top.loc[:, ~top.columns.isnull()]
    top.to_csv(args.out_top, sep="\t", index=False)
    print("Manhattan PNG:", args.out_png)
    print("Top hits:", args.out_top)

if __name__ == "__main__":
    main()

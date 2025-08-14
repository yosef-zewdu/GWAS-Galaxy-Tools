import argparse
import pandas as pd
import numpy as np
import statsmodels.api as sm
from statsmodels.tools.sm_exceptions import PerfectSeparationError

def compute_maf(col):
    nonmiss = col.dropna()
    n = nonmiss.shape[0]
    if n == 0:
        return np.nan
    alt_af = nonmiss.sum() / (2.0 * n)
    return min(alt_af, 1.0 - alt_af)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--genotypes", required=True)
    p.add_argument("--phenotypes_covariates", required=True)
    p.add_argument("--snp_annotation", default=None)
    p.add_argument("--covariates", default="")  # comma-separated covariate names
    p.add_argument("--min_maf", type=float, default=0.01)
    p.add_argument("--missing_code", default="NA")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    geno = pd.read_csv(args.genotypes, sep="\t", index_col=0)
    geno = geno.replace(args.missing_code, np.nan).astype(float)

    pheno = pd.read_csv(args.phenotypes_covariates, sep="\t")
    if 'sample_id' not in pheno.columns or 'phenotype' not in pheno.columns:
        raise ValueError("phenotypes_covariates.tsv must contain 'sample_id' and 'phenotype' columns")

    pheno = pheno.set_index('sample_id')
    cov_list = [c for c in args.covariates.split(",") if c.strip()]

    results = []
    for snp in geno.columns:
        g = geno[snp]
        # align samples
        df = pheno.join(g.rename('genotype'), how='inner')
        df = df.dropna(subset=['phenotype'])  # must have phenotype
        # For genotype, we'll keep rows where genotype is non-missing
        df = df[~df['genotype'].isna()].copy()
        n_non = df.shape[0]
        if n_non < 10:
            continue
        maf = compute_maf(df['genotype'])
        if pd.isna(maf) or maf < args.min_maf:
            continue

        # Prepare design matrix
        X = pd.DataFrame({'genotype': df['genotype']})
        for c in cov_list:
            if c not in df.columns:
                raise ValueError(f"Covariate '{c}' not found in phenotypes_covariates file")
            X[c] = df[c]
        X = sm.add_constant(X, has_constant='add')
        y = df['phenotype']
        try:
            model = sm.Logit(y, X).fit(disp=0, maxiter=100)
        except (PerfectSeparationError, np.linalg.LinAlgError, ValueError):
            # skip SNP if model fails
            continue
        # genotype coefficient
        if 'genotype' not in model.params:
            continue
        beta = model.params['genotype']
        se = model.bse['genotype']
        z = beta / se if se != 0 else np.nan
        pval = model.pvalues['genotype']
        row = {'snp_id': snp, 'beta': beta, 'SE': se, 'z': z, 'p': pval, 'MAF': maf, 'n': n_non}
        results.append(row)

    resdf = pd.DataFrame(results)
    if args.snp_annotation:
        ann = pd.read_csv(args.snp_annotation, sep="\t")
        resdf = resdf.merge(ann, on='snp_id', how='left')
        # reorder
        cols = ['snp_id','chrom','pos','beta','SE','z','p','MAF','n']
        resdf = resdf[cols]
    resdf.to_csv(args.out, sep="\t", index=False)
    print("Association results written to", args.out)

if __name__ == "__main__":
    main()

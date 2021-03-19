import pandas as pd
import numpy as np
import statsmodels.api as sm
import argparse
import scipy.stats 

merged_data = pd.read_csv(args.merged, sep='\t')

def residualize(group, covariates):
    """Regress expression on covariates and return residual."""
    df = group[["expression"] + covariates]
    result = sm.OLS(df["expression"], df[covariates]).fit().resid
    return(result)

def correlate(df, col_name):
    """Return p-value of correlation between residualized expression and genotype PCs."""
    return(scipy.stats.pearsonr(df["expression"], df[col_name])[1])

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--covariates_to_regress', nargs='*', default=None)
    parser.add_argument('--out')
    args = parser.parse_args()

    # Optionally residualize expression values on user-specified covariates
    if args.covariates_to_regress: 
        resid = merged_data.groupby("gene").apply(lambda group: 
            residualize(group, args.covariates_to_regress)).reset_index()
        merged_data = pd.merge(merged_data.drop(["gene", "expression"], axis=1), 
                               resid, left_index=True, right_on="level_1")
        merged_data = merged_data.set_index("level_1").rename({0: "expression"}, axis=1)

    # Create DataFrame of p-values for correlation between expression and all covariates
    cov_corr = pd.DataFrame()
    cov_columns = ["seq_center", "exam", "sex", "race", "global_ancestry", "local_ancestry"] + \
                  ["genotype_PC" + str(i) for i in range(10)] + \
                  ["expression_PC" + str(i) for i in range(10)]
    for col in cov_columns:
        tmp = merged_data.groupby("gene").apply(lambda group: correlate(group, col))
        cov_corr = pd.concat([cov_corr, tmp], axis=1).rename({0: "corr_" + col}, axis=1)
    cov_corr.to_csv(args.out, sep='\t')


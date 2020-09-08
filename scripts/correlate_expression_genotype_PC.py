import pandas as pd
import numpy as np
import statsmodels.api as sm
import argparse
import scipy.stats 

parser = argparse.ArgumentParser()
parser.add_argument('--merged')
parser.add_argument('--covariate_file')
parser.add_argument('--covariates_to_regress', nargs='*', default=None)
parser.add_argument('--out')

args = parser.parse_args()

merged_data = pd.read_csv(args.merged, sep='\t')

def residualize(group, covariates):
    df = group[["expression"] + covariates]
    result = sm.OLS(df["expression"], df[covariates]).fit().resid
    return(result)

if args.covariates_to_regress: # Residualize expression values on user-specified covariates
    resid = merged_data.groupby("gene").apply(lambda group: residualize(group, args.covariates_to_regress)).reset_index()
    resid = pd.merge(merged_data[["nwd_id"]], resid, left_index=True, right_on="level_1")
    resid = resid.set_index("level_1").rename({0: "expression"}, axis=1)
else:
    resid = merged_data[["nwd_id", "expression", "gene"]]

cov = pd.read_csv(args.covariate_file, sep='\t', index_col="nwd_id").drop("BroadUWPlate", axis=1)

# Return p-value of correlation between residualized expression and genotypic PCs
def correlate(resid_df, cov_df):
    merged = pd.merge(resid_df, cov_df, left_on="nwd_id", right_index=True)
    return(scipy.stats.pearsonr(merged["expression"], merged[cov_df.name])[1])

cov_corr = pd.DataFrame()
for col in cov.columns:
    tmp = resid.groupby("gene").apply(lambda group: correlate(group, cov[col]))
    cov_corr = pd.concat([cov_corr, tmp], axis=1).rename({0: "corr_" + col}, axis=1)

cov_corr.to_csv(args.out, sep='\t')

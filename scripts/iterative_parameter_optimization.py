import pandas as pd
import argparse
import numpy as np
import cvxpy as cp
import statsmodels.api as sm

def document_params(delta_file, betas_file, delta, betas):
    """Write current iteration of parameters to file."""
    with open(betas_file, 'w') as f:
        betas.to_csv(f, sep='\t')
    with open(delta_file, 'a') as f:
        f.write(str(delta))
        f.write('\n')

def calculate_residual(row, covariates):
    res = np.array(row["expression"] - row["genotype_Afr"] * row["beta_Afr"] - row["genotype_Eur"] * row["beta_Eur"])
    for cov in covariates:
        res = res - row["int_" + cov] * row[cov]
    res = res - (row["beta_Afr"] - row["beta_Eur"]) * row["genotype_Eur"] * row["race"]
    return(res)

def optimize_delta(df, covariates, unconstrained):
    """Optimize delta in likelihood model using quadratic programming solver."""
    b = np.array(df["expression"] - df["genotype_Afr"] * df["beta_Afr"] - df["genotype_Eur"] * df["beta_Eur"])
    for cov in covariates:
        b = b - df["int_" + cov] * df[cov]
    A = np.array((df["beta_Afr"] - df["beta_Eur"]) * df["genotype_Eur"] * df["race"])
    h = np.array([1, 0])
    P = np.matrix(np.dot(A, A))
    q = np.matrix(np.dot(A, -b))
    G = np.array([1, -1])
    delta = cp.Variable(1)
    if unconstrained:
        prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(delta, P) + q.T @ delta))
    else:
        prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(delta, P) + q.T @ delta), [G @ delta <= h])
    prob.solve()
    result = delta.value[0]
    return(result)

def optimize_betas(group, covariates):
    """Optimize betas/effects in likelihood model using ordinary least squares linear regression."""
    beta_Afr = group["genotype_Afr"] + group["delta"] * group["genotype_Eur"] * group["race"]
    beta_Eur = group["genotype_Eur"] - group["delta"] * group["genotype_Eur"] * group["race"]
    Y = group["expression"]
    df_dict = {"Y": Y, "beta_Afr": beta_Afr, "beta_Eur": beta_Eur}
    for cov in covariates:
        df_dict["int_" + cov] = group[cov]
    df = pd.DataFrame(df_dict)
    result = sm.OLS(df["Y"], df[["beta_Afr", "beta_Eur"] + ["int_" + cov for cov in covariates]]).fit().params
    return(result)

def update_params(df, covariates, betas=None, delta=None):
    """Update current parameters in provided dataframe."""
    if delta is not None:
        df["delta"] = delta
    elif betas is not None:
        col_drop_lst = ["beta_Afr", "beta_Eur"] + ["int_" + cov for cov in covariates]
        df = df.drop(columns=col_drop_lst)
        df = pd.merge(df, betas, left_on="gene", right_index=True)
    return(df)

def neg_control(group):
    """Exclude African-American ancestry-heterozygous individuals. Randomly assign a subset of 100 European-
       Americans to be a validation set; code as African-Americans in order to perform negative control."""
    val_idv = np.random.choice(group.loc[group.race == 0, "nwd_id"].values, size=100, replace=False)
    group_subset = group[-((group.race == 1) & (group.local_ancestry < 2))]
    group_subset["race"] = group_subset.apply(lambda row: 1 if row.nwd_id in val_idv else row.race, axis=1, result_type='expand')
    return(group_subset)

def remove_asc(group, partition_matrix):
    """Remove European ascertainment individuals from dataset for optimizing likelihood model. These 
       individuals were used to select the genes/variants we are optimizing the model on, so inclusion
       of these individuals will result in inflated European effect sizes."""
    gene = group.name
    asc_idv = partition_matrix.loc[gene, partition_matrix.loc[gene] == 1].index
    return(group[-group.nwd_id.isin(asc_idv)])

def bootstrap_data(all_df):
    """Generate one random sample (with replacement) from the input dataframe. Size of sample is 
       specified by number of genes in input dataframe."""
    def decrement(d): # No return statement because dictionaries are modified in-place
        keys_to_drop = []
        for key, val in d.items():
            val = val - 1
            if val == 0:
                keys_to_drop.append(key)
            else:
                d[key] = val
        for key in keys_to_drop:
            del d[key]

    gene_list = all_df["gene"].unique()
    n_genes = len(gene_list)
    bootstrap_list = np.random.choice(gene_list, n_genes, replace=True)
    bootstrap_dict = {}
    for gene in bootstrap_list:
        bootstrap_dict[gene] = bootstrap_dict.get(gene, 0) + 1

    bootstrap_df = pd.DataFrame()
    while bootstrap_dict:
        tmp_df = all_df[all_df["gene"].isin(bootstrap_dict.keys())]
        bootstrap_df = pd.concat([bootstrap_df, tmp_df])
        decrement(bootstrap_dict)

    return(bootstrap_df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--max_iter')
    parser.add_argument('--delta_out')
    parser.add_argument('--betas_out')
    parser.add_argument('--group')
    parser.add_argument('--partition')
    parser.add_argument('--covariates', nargs='*')
    parser.add_argument('--delta', default=None)
    parser.add_argument('--residuals', default=None)
    parser.add_argument('--bootstrap', action='store_true')
    parser.add_argument('--unconstrained', action='store_true')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')
    merged_data = merged_data[args.covariates + ["genotype_Afr", "genotype_Eur", "expression", "nwd_id", "gene", "race"]]
    
    # Initialize empty parameter columns
    merged_data["beta_Afr"] = None
    merged_data["beta_Eur"] = None
    for cov in args.covariates:
        merged_data["int_" + cov] = None

    # Remove ascertainment individuals to avoid biasing parameter estimation
    partition_matrix = pd.read_csv(args.partition, sep='\t', index_col=0)
    merged_data = merged_data.groupby("gene").apply(lambda grp: remove_asc(grp, partition_matrix)).reset_index(drop=True)

    # Deal with command line arguments/options
    if args.group == "control":
        merged_data = merged_data.groupby("gene").apply(neg_control).reset_index(drop=True)
    if args.bootstrap:
        merged_data = bootstrap_data(merged_data)

    if args.delta is not None: # Optimize betas for a user-specified delta + write to file
        curr_delta = args.delta
        merged_data = update_params(merged_data, args.covariates, delta=curr_delta)
        curr_betas = merged_data.groupby("gene").apply(lambda group: optimize_betas(group, args.covariates))
        merged_data = update_params(merged_data, args.covariates, betas=curr_betas)
        document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)
    else: # Iteratively optimize delta + betas
        prev_delta = -1 # Ensures we do 1+ iterations
        for i in range(int(args.max_iter)):
            if i == 0:
                curr_delta = np.random.uniform(0, 1)
            else:
                curr_delta = optimize_delta(merged_data, args.covariates, args.unconstrained)
            if (abs(curr_delta - prev_delta) < .0001):
                break
            merged_data = update_params(merged_data, args.covariates, delta=curr_delta)
            prev_delta = curr_delta
            curr_betas = merged_data.groupby("gene").apply(lambda group: optimize_betas(group, args.covariates))
            merged_data = update_params(merged_data, args.covariates, betas=curr_betas)
            document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)

    if args.residuals is not None: # Write residuals to file
        fitted_residuals = merged_data.apply(lambda row: calculate_residual(row, args.covariates), axis=1)
        fitted_residuals.to_csv(args.residuals, index=False)


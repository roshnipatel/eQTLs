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

def optimize_delta(df, covariates):
    """Optimize delta in likelihood model using ordinary least squares 
       linear regression."""
    Y = np.array(df["expression"] - df["genotype_Afr"] * df["beta_Afr"] - 
                 df["genotype_Eur"] * df["beta_Eur"])
    for cov in covariates:
        Y = Y - df["int_" + cov] * df[cov]
    X = np.array((df["beta_Afr"] - df["beta_Eur"]) * 
                 df["genotype_Eur"] * df["race"])
    result = sm.OLS(Y, X).fit().params[0]
    return(result)

def optimize_betas(group, covariates):
    """Optimize betas/effects in likelihood model using ordinary least squares 
       linear regression."""
    beta_Afr = group["genotype_Afr"] + group["delta"] * group["genotype_Eur"] * group["race"]
    beta_Eur = group["genotype_Eur"] - group["delta"] * group["genotype_Eur"] * group["race"]
    Y = group["expression"]
    df_dict = {"Y": Y, "beta_Afr": beta_Afr, "beta_Eur": beta_Eur}
    for cov in covariates:
        df_dict["int_" + cov] = group[cov]
    df = pd.DataFrame(df_dict)
    result = sm.OLS(df["Y"], df[["beta_Afr", "beta_Eur"] + 
                    ["int_" + cov for cov in covariates]]).fit().params
    return(result)

def fit_no_beta(group, covariates):
    """Optimize likelihood model without fitting SNP beta/effect using 
       ordinary least squares linear regression."""
    Y = group["expression"]
    df_dict = {"Y": Y}
    for cov in covariates:
        df_dict["int_" + cov] = group[cov]
    df = pd.DataFrame(df_dict)
    if ("race_Afr" in covariates) and ("race_Eur" in covariates):
        result = sm.OLS(df["Y"], 
                        df[["int_" + cov for cov in covariates]]).fit().params
    else:
        df["intercept"] = 1
        result = sm.OLS(df["Y"], df[["intercept"] + 
                        ["int_" + cov for cov in covariates]]).fit().params
    return(result)

def fit_single_beta(group, covariates):
    """Optimize single beta in likelihood model using ordinary least squares 
       linear regression."""
    beta = group["genotype_Afr"] + group["genotype_Eur"]
    Y = group["expression"]
    df_dict = {"Y": Y, "beta": beta}
    for cov in covariates:
        df_dict["int_" + cov] = group[cov]
    df = pd.DataFrame(df_dict)
    if ("race_Afr" in covariates) and ("race_Eur" in covariates):
        result = sm.OLS(df["Y"], df[["beta"] + 
                        ["int_" + cov for cov in covariates]]).fit().params
    else:
        df["intercept"] = 1
        result = sm.OLS(df["Y"], df[["beta", "intercept"] + 
                        ["int_" + cov for cov in covariates]]).fit().params
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
    """Exclude African-American ancestry-heterozygous individuals. Randomly 
       assign a subset of 100 European-Americans to be a validation set; 
       code as African-Americans in order to perform negative control."""
    val_idv = np.random.choice(group.loc[group.race == 0, "nwd_id"].values, 
                               size=100, replace=False)
    group_subset = group[-((group.race == 1) & (group.local_ancestry < 2))]
    group_subset["race"] = group_subset.apply(lambda row: 1 if row.nwd_id in val_idv 
                                              else row.race, axis=1, result_type='expand')
    return(group_subset)

def remove_ind(group, partition_matrix):
    """Remove European ascertainment individuals from dataset for optimizing 
       likelihood model. These individuals were used to select the genes/variants 
       we are optimizing the model on, so inclusion of these individuals will 
       result in inflated European effect sizes."""
    gene = group.name
    asc_idv = partition_matrix.loc[gene, partition_matrix.loc[gene] == 1].index
    return(group[-group.nwd_id.isin(asc_idv)])

def bootstrap(df):


def bootstrap_over_snps(all_df):
    """Generate one random sample (with replacement) from the input dataframe. 
       Size of sample is specified by number of genes in input dataframe."""
    def decrement(d): 
        """For all keys in input dictionary, decrement values by 1. If value equals
           0, remove key from dictionary. Function has no return statement 
           because dictionaries are modified in-place."""
        keys_to_drop = []
        for key, val in d.items():
            val = val - 1
            if val == 0:
                keys_to_drop.append(key)
            else:
                d[key] = val
        for key in keys_to_drop:
            del d[key]

    # Randomly sample over genes with replacement and document number of samples 
    # for each gene in dictionary
    gene_list = all_df["gene"].unique()
    n_genes = len(gene_list)
    bootstrap_list = np.random.choice(gene_list, n_genes, replace=True)
    bootstrap_dict = {}
    for gene in bootstrap_list:
        bootstrap_dict[gene] = bootstrap_dict.get(gene, 0) + 1

    # Subset full dataset for all genes in bootstrap dictionary. Then decrement
    # values for all genes in dictionary, removing gene if value reaches 0.
    # Repeat until dictionary is empty.
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
    parser.add_argument('--mode')
    parser.add_argument('--ascertainment')
    parser.add_argument('--idv_bin', default=None)
    parser.add_argument('--gene_bin', default=None)
    parser.add_argument('--validation', default=None)
    parser.add_argument('--covariates', nargs='*')
    parser.add_argument('--delta', default=None)
    parser.add_argument('--bootstrap', action='store_true')
    parser.add_argument('--jackknife')
    parser.add_argument('--single_beta', action='store_true')
    parser.add_argument('--no_beta', action='store_true')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')
    col_names = args.covariates + ["genotype_Afr", "genotype_Eur", "expression", 
                                   "nwd_id", "gene", "race"]
    merged_data = merged_data[col_names]

    # Optionally filter for the individuals and/or genes in the provided bin file
    if args.idv_bin:
        idv_bin = pd.read_csv(args.idv_bin, squeeze=True)
        merged_data = merged_data.loc[merged_data.nwd_id.isin(idv_bin)] 
    if args.gene_bin:
        gene_bin = pd.read_csv(args.gene_bin, squeeze=True)
        merged_data["trimmed_gene"] = merged_data.apply(lambda row: row.gene[:15], axis=1)
        merged_data = merged_data.loc[merged_data.trimmed_gene.isin(gene_bin)]
        merged_data = merged_data.drop(columns=["trimmed_gene"])

    # If ascertainment matrix provided, remove those individuals to avoid biasing 
    # parameter estimation (provided for all executions of script except on 
    # simulated data.
    if args.ascertainment:
        ascertainment_matrix = pd.read_csv(args.ascertainment, sep='\t', index_col=0)
        merged_data = merged_data.groupby("gene").apply(lambda grp: 
            remove_ind(grp, ascertainment_matrix)).reset_index(drop=True)
    print("removed asc")

    # If estimating parameters for negative control, massage data accordingly
    if args.group == "control":
        merged_data = merged_data.groupby("gene").apply(neg_control).reset_index(drop=True)

    # Optionally bootstrap over data, randomly sampling genes with replacement
    if args.bootstrap:
        merged_data = bootstrap(merged_data)

    # Optionally perform jackknife, dropping the user-specified gene from the data
    if args.jackknife:
        merged_data = merged_data.loc[merged_data.gene != args.jackknife, :]

    # If validation matrix provided, remove those individuals before estimating parameters
    if args.validation:
        validation_matrix = pd.read_csv(args.validation, sep='\t', index_col=0)
        merged_data = merged_data.groupby("gene").apply(lambda grp: 
            remove_ind(grp, validation_matrix)).reset_index(drop=True)
    print("removed val")
    
    # Initialize empty parameter columns
    merged_data["beta_Afr"] = None
    merged_data["beta_Eur"] = None
    for cov in args.covariates:
        merged_data["int_" + cov] = None

    # Optimize parameters
    if args.delta is not None: # Fit user-specified delta
        curr_delta = float(args.delta)
        merged_data = update_params(merged_data, args.covariates, delta=curr_delta)
        curr_betas = merged_data.groupby("gene").apply(lambda group: optimize_betas(group, args.covariates))
        merged_data = update_params(merged_data, args.covariates, betas=curr_betas)
        document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)
    elif args.mode != "fit_delta": # Fit model without delta
        curr_delta = 0
        merged_data = update_params(merged_data, args.covariates, delta=curr_delta)
        if args.mode == "fit_no_beta":
            curr_betas = merged_data.groupby("gene").apply(lambda group: fit_no_beta(group, args.covariates))
        elif args.mode == "fit_single_beta":
            curr_betas = merged_data.groupby("gene").apply(lambda group: fit_single_beta(group, args.covariates))
        elif args.mode == "no_interaction":
            curr_betas = merged_data.groupby("gene").apply(lambda group: optimize_betas(group, args.covariates))
        merged_data = update_params(merged_data, args.covariates, betas=curr_betas)
        document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)
    elif args.mode == "fit_delta": # Fit model with delta, iteratively optimizing delta + betas
        print("fitting delta")
        prev_delta = -1 # Ensures we do 1+ iterations
        for i in range(int(args.max_iter)):
            if i == 0:
                curr_delta = np.random.uniform(0, 1)
            else:
                curr_delta = optimize_delta(merged_data, args.covariates)
            if (abs(curr_delta - prev_delta) < .0001):
                break
            merged_data = update_params(merged_data, args.covariates, delta=curr_delta)
            prev_delta = curr_delta
            curr_betas = merged_data.groupby("gene").apply(lambda group: optimize_betas(group, args.covariates))
            merged_data = update_params(merged_data, args.covariates, betas=curr_betas)
            document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)


import pandas as pd
import numpy as np
import argparse
from iterative_parameter_optimization import remove_ind

def compute_var(curr_betas, delta, exp_df, model_terms):
    """For a single gene, compute the proportion of 
       expression variance unexplained by the model."""
    gene = curr_betas.gene
    curr_gene = exp_df.loc[exp_df.gene == gene]
    # Sum up model terms
    mod = np.zeros(len(curr_gene))
    for term in model_terms:
        mod += curr_betas[term] * curr_gene.loc[:,term]
    if delta != 0:
        mod += delta * (curr_betas.genotype_Afr - curr_betas.genotype_Eur) * \
                curr_gene.genotype_Eur * curr_gene.race
    var_exp = np.var(curr_gene.loc[:,"expression"])
    resid = curr_gene.loc[:,"expression"] - mod
    var_resid = np.var(resid)
    # Return proportion residual variance
    return(var_resid / var_exp)

def maf_filtered_genes(merged, maf_threshold):
    def filter_helper(gene_df, maf_df, threshold):
        filt_genes = gene_df.loc[(maf_df > threshold) & (maf_df < 1 - threshold)].reset_index().gene
        return(filt_genes)
    sum_df = merged.groupby(["gene", "ID"]).sum()
    afr_maf = sum_df.genotype_Afr / sum_df.local_ancestry 
    eur_maf = sum_df.genotype_Eur / ((sum_df.race_Afr + sum_df.race_Eur) * 2 - sum_df.local_ancestry)
    afr_genes = filter_helper(sum_df, afr_maf, maf_threshold)
    eur_genes = filter_helper(sum_df, eur_maf, maf_threshold)
    return(pd.merge(afr_genes, eur_genes))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--delta')
    parser.add_argument('--betas')
    parser.add_argument('--ind')
    parser.add_argument('--validation')
    parser.add_argument('--model_terms', nargs='+')
    parser.add_argument('--mode')
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')

    # Remove non-validation individuals to avoid biasing variance estimates
    validation_matrix = pd.read_csv(args.validation, sep='\t', index_col=0)
    remove_matrix = 1 - validation_matrix
    merged_data = merged_data.groupby("gene").apply(lambda grp: 
                  remove_ind(grp, remove_matrix)).reset_index(drop=True)

    # Perform MAF filtering within smaller set of validation individuals
    filtered_genes = maf_filtered_genes(merged_data, 0.05)

    # Read in fitted delta
    with open(args.delta, 'r') as f:
        delta = float(f.readlines()[-1].strip())

    # Read in fitted betas
    betas = pd.read_csv(args.betas, sep='\t', index_col=0)
    betas = pd.merge(betas, filtered_genes, left_index=True, right_on="gene")

    # Rename columns of betas to match columns of merged_data
    betas_cols = betas.columns
    rename_dict = {}
    for col in betas_cols:
        if col[:4] == "int_":
            rename_dict[col] = col[4:]
        elif col[:5] == "beta_":
            rename_dict[col] = "genotype_" + col[5:]
        else:
            rename_dict[col] = col
    betas = betas.rename(columns=rename_dict)
    
    # If using a single beta b (i.e. not ancestry-specific betas bA and bE)
    # then set coefficients of gA and gE equal to the same beta, b
    if "beta" in rename_dict:
        betas["genotype_Afr"] = betas.loc[:,"beta"]
        betas["genotype_Eur"] = betas.loc[:,"beta"]

    # If using model with single intercept c (i.e. not race-specific 
    # intercepts cA and cE) then set coefficients of rA and rE equal to the
    # same intercept, c
    if "intercept" in rename_dict:
        betas["race_Afr"] = betas.loc[:,"intercept"]
        betas["race_Eur"] = betas.loc[:,"intercept"]

    # If fitting genotype betas, add these terms to the set model_terms 
    model_terms = set(args.model_terms)
    model_terms.add("race_Afr")
    model_terms.add("race_Eur")
    if args.mode != "fit_no_beta":
        model_terms.add("genotype_Afr") 
        model_terms.add("genotype_Eur") 

    # Compute proportion variance explained by model in validation set
    betas.loc[:,"prop_var"] = betas.apply(lambda row: 
        compute_var(row, delta, merged_data, model_terms), axis=1)
    betas.to_csv(args.out, index=False, columns=["gene", "prop_var"], sep='\t')

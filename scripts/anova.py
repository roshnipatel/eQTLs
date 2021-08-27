import pandas as pd
import numpy as np
import argparse
from iterative_parameter_optimization import remove_ind

def compute_var(curr_betas, delta, exp_df, model_terms, resid_terms, ind):
    """For a single gene, compute the proportion of (optionally residualized)
       expression variance explained by the model."""
    gene = curr_betas.gene
    
    # Specify which individuals are used to calculate variance explained
    if ind == "all": # Use all individuals
        curr_gene = exp_df.loc[exp_df.gene == gene]
    elif ind == "AAEur": # Use only AA with Eur ancestry at locus
        curr_gene = exp_df.loc[(exp_df.gene == gene) & ((exp_df.race == 1) 
                               & (exp_df.local_ancestry < 2))]
    elif ind == "AAAfr": # Use only AA with Afr ancestry at locus
        curr_gene = exp_df.loc[(exp_df.gene == gene) & ((exp_df.race == 1) 
                               & (exp_df.local_ancestry == 2))]
    elif ind == "EA": # Use only EA
        curr_gene = exp_df.loc[(exp_df.gene == gene) & (exp_df.race == 0)]
    elif ind == "AA": # Use only AA
        curr_gene = exp_df.loc[(exp_df.gene == gene) & (exp_df.race == 1)]

    # Sum up model terms
    mod = np.zeros(len(curr_gene))
    for term in model_terms:
        mod = curr_betas[term] * curr_gene.loc[:,term]
    if delta != 0:
        mod += delta * (curr_betas.genotype_Afr - curr_betas.genotype_Eur) * \
                curr_gene.genotype_Eur * curr_gene.race
   
    # Sum up terms to residualize
    resid = np.zeros(len(curr_gene))
    for term in resid_terms:
        resid = curr_betas[term] * curr_gene.loc[:,term]

    var_exp = np.var(curr_gene.loc[:,"expression"] - resid)
    var_mod = np.var(mod)

    # Return proportion variance explained by model 
    return(var_mod / var_exp)

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
    parser.add_argument('--residualized_terms', nargs='*')
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
    if "beta" in rename_dict:
        betas["genotype_Afr"] = betas.loc[:,"beta"]
        betas["genotype_Eur"] = betas.loc[:,"beta"]

    # Parse model terms and terms to residualize, inferring additional
    # model terms from analysis mode
    model_terms = args.model_terms
    if len(model_terms) == 2 and args.mode == "fit_no_beta": 
        model_terms.append("intercept")
        merged_data["intercept"] = 1
    elif args.mode != "fit_no_beta":
        model_terms.extend(["genotype_Afr", "genotype_Eur"]) 
    resid_terms = args.residualized_terms

    # Compute proportion variance explained by model in validation set
    betas["prop_var"] = betas.apply(lambda row: 
        compute_var(row, delta, merged_data, model_terms, resid_terms, args.ind), axis=1)
    betas.to_csv(args.out, index=True, columns=["prop_var"], sep='\t')

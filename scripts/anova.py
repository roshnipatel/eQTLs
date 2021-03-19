import pandas as pd
import numpy as np
import argparse
from iterative_parameter_optimization import remove_ind

def compute_var(exp_df, delta, coeff, curr_betas, ind):
    """For a single gene, compute the proportion of expression variance 
       explained by the model."""
    gene = curr_betas.name
    
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

    # Calculate variance of gene expression
    var_exp = np.var(curr_gene.loc[:,"expression"])

    # Calculate variance of predicted expression
    pred = np.zeros(len(curr_gene))
    for c in coeff:
        if c == "intercept":
            pred += curr_betas[c]
        elif c[:3] == "int":
            pred += curr_betas[c] * curr_gene.loc[:,c[4:]] 
        elif c == "beta":
            pred += curr_betas[c] * (curr_gene.loc[:,"genotype_Afr"] + 
                    curr_gene.loc[:,"genotype_Eur"]) 
        else:
            pred += curr_betas[c] * curr_gene.loc[:,"genotype_" + c[-3:]] 
    if delta != 0:
        pred += delta * (curr_betas.beta_Afr - curr_betas.beta_Eur) * \
                curr_gene.genotype_Eur * curr_gene.race
    var_pred = np.var(pred)

    # Return proportion variance explained by model 
    return(var_pred / var_exp)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--delta')
    parser.add_argument('--betas')
    parser.add_argument('--ind')
    parser.add_argument('--validation')
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')

    # Remove non-validation individuals to avoid biasing variance estimates
    validation_matrix = pd.read_csv(args.validation, sep='\t', index_col=0)
    remove_matrix = 1 - validation_matrix
    merged_data = merged_data.groupby("gene").apply(lambda grp: 
                  remove_ind(grp, remove_matrix)).reset_index(drop=True)

    # Read in fitted model coefficients (betas for each gene and delta)
    betas = pd.read_csv(args.betas, sep='\t', index_col=0)
    model_coeff = betas.columns
    with open(args.delta, 'r') as f:
        delta = float(f.readlines()[-1].strip())

    # Compute proportion variance explained by model in validation set
    betas["prop_var"] = betas.apply(lambda row: 
           compute_var(merged_data, delta, model_coeff, row, args.ind), axis=1)
    betas.to_csv(args.out, index=True, columns=["prop_var"], sep='\t')

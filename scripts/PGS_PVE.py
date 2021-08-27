import pandas as pd
import numpy as np
import argparse
from iterative_parameter_optimization import remove_ind
from anova import compute_var, maf_filtered_genes

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--delta')
    parser.add_argument('--betas')
    parser.add_argument('--ind')
    parser.add_argument('--validation')
    parser.add_argument('--residualized_terms', nargs='*')
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
    print(betas_cols)
    rename_dict = {}
    for col in betas_cols:
        if col[:4] == "int_":
            rename_dict[col] = col[4:]
        elif col[:5] == "beta_":
            rename_dict[col] = "genotype_" + col[5:]
    print(rename_dict)
    betas = betas.rename(columns=rename_dict)

    # Create df of "PGS" betas where we use European beta (denoted by
    # genotype_Eur column) to approximate SNP effect in African ancestry
    PGS_betas = betas.copy()
    PGS_betas["genotype_Afr"] = PGS_betas["genotype_Eur"]

    model_terms = ["genotype_Eur", "genotype_Afr"]
    resid_terms = args.residualized_terms
    betas["PGS_prop_var"] = PGS_betas.apply(lambda row:
        compute_var(row, 0, merged_data, model_terms, resid_terms, args.ind), axis=1)
    betas["full_model_prop_var"] = betas.apply(lambda row:
        compute_var(row, delta, merged_data, model_terms, resid_terms, args.ind), axis=1)
    betas.to_csv(args.out, index=True, 
                 columns=["PGS_prop_var", "full_model_prop_var"], sep='\t')

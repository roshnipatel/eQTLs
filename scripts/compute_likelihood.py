import pandas as pd
import argparse
import numpy as np
from iterative_parameter_optimization import update_params
from iterative_parameter_optimization import optimize_betas
from iterative_parameter_optimization import drop_asc, neg_control

def likelihood(df, delta, covariates):
    """Returns log-likelihood of data based on the given value of delta."""
    def rowwise_likelihood(row, delta, covariates):
        exp = row.expression
        additive_effects = row.genotype_Afr * row.beta_Afr + row.genotype_Eur * row.beta_Eur
        interaction = delta * (row.beta_Afr - row.beta_Eur) * row.genotype_Eur * row.race
        coeff = 0
        for cov in args.covariates:
            coeff += row["int_" + cov] * row[cov]
        row_likelihood = -(exp - additive_effects - coeff - interaction) ** 2
        return(row_likelihood)
    df_likelihood = df.apply(lambda x: rowwise_likelihood(x, delta, covariates), axis=1)
    return(df_likelihood.sum())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--option')
    parser.add_argument('--betas', default=None)
    parser.add_argument('--delta', default=None)
    parser.add_argument('--covariates', nargs='+')
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')

    # Deal with command line arguments/options
    if args.option == "neg_control":
        merged_data = merged_data.groupby("gene").apply(neg_control).reset_index(drop=True)
    if args.option == "drop_asc":
        merged_data = merged_data.groupby("gene").apply(drop_asc).reset_index(drop=True)
    
    # Initialize empty parameter columns
    merged_data["beta_Afr"] = None
    merged_data["beta_Eur"] = None
    for cov in args.covariates:
        merged_data["int_" + cov] = None

    if args.delta is not None: # Compute likelihood of data only for user-specified value of delta
        delta_list = [float(args.delta)]
    else: # Compute likelihood of data for 100 values of delta over a uniform grid on [0, 1]
        delta_list = [i / 100 for i in range(0, 101)]

    if args.betas is not None: # Compute likelihood conditional on user-specified betas
        betas = pd.read_csv(args.betas, sep='\t')
        betas = betas.set_index("gene")
        merged_data = update_params(merged_data, args.covariates, betas=betas)

    likelihood_list = []
    for d in delta_list:
        print("processing delta = {0}".format(str(d)))
        if args.betas is None: # Optimize betas for each value of delta, and then compute likelihood
            merged_data = update_params(merged_data, args.covariates, delta=d)
            curr_betas = merged_data.groupby("gene").apply(lambda group: optimize_betas(group, args.covariates))
            merged_data = update_params(merged_data, args.covariates, betas=curr_betas)
        likelihood_list.append(likelihood(merged_data, d, args.covariates))
    res = pd.DataFrame({"Delta": delta_list, "Likelihood": likelihood_list})
    res.to_csv(args.out, sep='\t', index=False)

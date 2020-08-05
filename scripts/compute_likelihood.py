import pandas as pd
import argparse
import numpy as np
from iterative_parameter_optimization import update_params
from iterative_parameter_optimization import drop_asc

def likelihood(df, delta, covariates):
    """Returns log-likelihood of data based on the given value of delta."""
    def rowwise_likelihood(row, delta, covariates):
        exp = row.Expression
        additive_effects = row.Genotype_Afr * row.Beta_Afr + row.Genotype_Eur * row.Beta_Eur
        interaction = delta * (row.Beta_Afr - row.Beta_Eur) * row.Genotype_Eur * row.Race_AA
        coeff = 0
        for cov in args.covariates:
            coeff += row["Int_" + cov] * row[cov]
        row_likelihood = -(exp - additive_effects - coeff - interaction) ** 2
        return(row_likelihood)
    df_likelihood = df.apply(lambda x: rowwise_likelihood(x, delta, covariates), axis=1)
    return(df_likelihood.sum())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--drop_asc', action='store_true')
    parser.add_argument('--betas', default=None)
    parser.add_argument('--delta', default=None)
    parser.add_argument('--covariates', nargs='+')
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')
    
    # Initialize empty parameter columns
    merged_data["Beta_Afr"] = None
    merged_data["Beta_Eur"] = None
    for cov in args.covariates:
        merged_data["Int_" + cov] = None

    # Remove ascertainment individuals to avoid artificial inflation of likelihood
    if args.drop_asc: 
        merged_data = merged_data.groupby("Gene").apply(drop_asc).reset_index(drop=True)

    if args.betas is not None:
        betas = pd.read_csv(args.betas, sep='\t')
        betas = betas.set_index("Gene")
        merged_data = update_params(merged_data, args.covariates, betas=betas)

    if args.delta is not None: # Compute likelihood of data for user-specified value of delta
        d_likelihood = likelihood(merged_data, float(args.delta), args.covariates)
        with open(args.out, 'w') as f:
            f.write(str(d_likelihood))
            f.write('\n')
    else: # Compute likelihood of data for 100 values of delta over a uniform grid on [0, 1]
        delta_list = [i / 100 for i in range(0, 101)]
        likelihood_list = []
        for d in delta_list:
            print("processing delta = {0}".format(str(d)))
            likelihood_list.append(likelihood(merged_data, d))
        res = pd.DataFrame({"Delta": delta_list, "Likelihood": likelihood_list})
        res.to_csv(args.out, sep='\t', index=False)


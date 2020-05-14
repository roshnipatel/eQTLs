import pandas as pd
import argparse
import numpy as np
from iterative_parameter_optimization import update_params

def likelihood(df, delta):
    def rowwise_likelihood(row, delta):
        exp = row.Expression
        additive_effects = row.Genotype_Afr * row.Curr_Effect_Afr + row.Genotype_Eur * row.Curr_Effect_Eur
        interaction = delta * (row.Curr_Effect_Afr - row.Curr_Effect_Eur) * row.Genotype_Eur * row.Race_AA
        row_likelihood = (exp - additive_effects - interaction) ** 2
        return(row_likelihood)
    df_likelihood = df.apply(lambda x: rowwise_likelihood(x, delta), axis=1)
    return(df_likelihood.sum())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--betas', default=None)
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')
    merged_data[["Curr_Effect_Afr", "Curr_Effect_Eur"]] = merged_data[["effect_Afr", "effect_Eur"]]
    if args.betas is not None:
        betas = pd.read_csv(args.betas, sep='\t')
        betas = betas.set_index("Gene")
        merged_data = update_params(merged_data, betas=betas)

    # Compute likelihood of data for a range of delta values, assuming fixed betas
    delta_list = [i / 100 for i in range(0, 101)]
    likelihood_list = []
    for d in delta_list:
        print("processing delta = {0}".format(str(d)))
        likelihood_list.append(likelihood(merged_data, d))
    res = pd.DataFrame({"Delta": delta_list, "Likelihood": likelihood_list})
    res.to_csv(args.out, sep='\t', index=False)


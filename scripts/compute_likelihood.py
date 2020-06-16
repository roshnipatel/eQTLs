import pandas as pd
import argparse
import numpy as np
from iterative_parameter_optimization import update_params
from iterative_parameter_optimization import drop_asc

def likelihood(df, delta):
    ### Note that this actually gives you the NEGATIVE log likelihood and I just haven't changed so as not to confuse myself
    def rowwise_likelihood(row, delta):
        exp = row.Expression
        additive_effects = row.Genotype_Afr * row.Curr_Effect_Afr + row.Genotype_Eur * row.Curr_Effect_Eur + row.Race_AA * row.Curr_Int_Afr + (1 - row.Race_AA) * row.Curr_Int_Eur
        interaction = delta * (row.Curr_Effect_Afr - row.Curr_Effect_Eur) * row.Genotype_Eur * row.Race_AA
        row_likelihood = (exp - additive_effects - interaction) ** 2
        return(row_likelihood)
    df_likelihood = df.apply(lambda x: rowwise_likelihood(x, delta), axis=1)
    return(df_likelihood.sum())

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--drop_asc', action='store_true')
    parser.add_argument('--betas', default=None)
    parser.add_argument('--delta', default=None)
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = pd.read_csv(args.merged, sep='\t')
    merged_data[["Curr_Effect_Afr", "Curr_Effect_Eur"]] = merged_data[["effect_Afr", "effect_Eur"]]
    merged_data[["Curr_Int_Afr", "Curr_Int_Eur"]] = merged_data[["intercept_Afr", "intercept_Eur"]]
    if args.drop_asc:
        merged_data = merged_data.groupby("Gene").apply(drop_asc).reset_index(drop=True)

    if args.betas is not None:
        betas = pd.read_csv(args.betas, sep='\t')
        betas = betas.set_index("Gene")
        merged_data = update_params(merged_data, betas=betas)

    if args.delta is not None:
        d_likelihood = likelihood(merged_data, float(args.delta))
        with open(args.out, 'w') as f:
            f.write(str(d_likelihood))
            f.write('\n')
    else:
        delta_list = [i / 100 for i in range(0, 101)]
        likelihood_list = []
        for d in delta_list:
            print("processing delta = {0}".format(str(d)))
            likelihood_list.append(likelihood(merged_data, d))
        res = pd.DataFrame({"Delta": delta_list, "Likelihood": likelihood_list})
        res.to_csv(args.out, sep='\t', index=False)


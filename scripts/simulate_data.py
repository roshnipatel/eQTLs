import pandas as pd
import numpy as np
import argparse
import scipy.stats

# Specify parameters for the beta distribution used to sample global ancestry
SHAPE_1 = 7.9
SHAPE_2 = 2.1

def simulate_ganc(row):
    """Simulate global ancestry from beta distribution specified by shape
       parameters at top of script. Only simulates global African ancestry
       greater than 0.5, since our analysis pipeline filters individuals
       as such."""
    ganc = 0
    if row[0] == 1:
        while ganc < 0.5:
            ganc = np.random.beta(SHAPE_1, SHAPE_2)
    return(ganc)

def simulate_lanc(row, n):
    """Simulate local ancestry as binomial conditional on race/population."""
    lanc = np.random.binomial(2, row[0], size=n)
    return(lanc)

def simulate_genotypes(row, SNP_df):
    """Simulate genotypes as binomial conditional on local ancestry."""
    geno = []
    for lanc in row:
        geno_Afr = np.random.binomial(lanc, SNP_df.loc[row.name, "MAF_Afr"])
        geno_Eur = np.random.binomial(2 - lanc, SNP_df.loc[row.name, "MAF_Eur"])
        geno.append(str(geno_Eur) + "-" + str(geno_Afr))
    return(geno)

def simulate_phenotypes(row, SNP_df, race_df, params):
    """Simulate phenotypes by multiplying genotypes by betas and adding random, 
       normally-distributed noise."""
    pheno = []
    beta_Afr = SNP_df.loc[row.name, "Beta_Afr"]
    beta_Eur = SNP_df.loc[row.name, "Beta_Eur"]
    for i in range(len(row)):
        geno = row.values[i]
        ind = row.index[i]
        ind_race = race_df.loc[ind, "race"]
        geno_Eur = int(geno[0])
        geno_Afr = int(geno[2])
        pheno_Afr = geno_Afr * beta_Afr
        pheno_Eur = geno_Eur * beta_Eur
        if ind_race == 1:
            noise = np.random.normal(0, np.sqrt(params["var_error_afr"]))
        else:
            noise = np.random.normal(0, np.sqrt(params["var_error_eur"]))
        interaction = params["delta"] * (beta_Afr - beta_Eur) * geno_Eur * ind_race
        pheno.append(pheno_Afr + pheno_Eur + noise + interaction)
    return(pheno)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--sim_betas', help="table of simulated betas", default=None)
    parser.add_argument('--error_afr', help="specified error term", default=None)
    parser.add_argument('--error_eur', help="specified error term", default=None)
    parser.add_argument('--delta', help="specified delta term", default=None)
    parser.add_argument('--out', help="output")
    args = parser.parse_args()

    # Parse command-line arguments and specify parameters for simulation
    params = {"n_idv": 200}
    if args.error_afr:
        params["var_error_afr"] = float(args.error_afr)
    else:
        params["var_error_afr"] = 0.15
    if args.error_eur:
        params["var_error_eur"] = float(args.error_eur)
    else:
        params["var_error_eur"] = 0.15
    if args.delta:
        params["delta"] = float(args.delta)
    else:
        params["delta"] = 0

    sig_SNPs = pd.read_csv(args.sim_betas, sep='\t')
    n_sig_SNPs = sig_SNPs.shape[0]

    # Simulate local ancestry, genotypes, and phenotypes for each of the 
    # significant SNPs
    race = pd.DataFrame([anc for anc in [0, 1] for _ in range(params["n_idv"])], 
                        columns=["race"])
    ganc = pd.DataFrame(race.apply(simulate_ganc, axis=1))
    lanc = ganc.apply(lambda x: simulate_lanc(x, n_sig_SNPs), axis=1, 
                      result_type='expand').T
    geno = lanc.apply(lambda x: simulate_genotypes(x, sig_SNPs), axis=1, 
                      result_type='expand')
    pheno = geno.apply(lambda x: simulate_phenotypes(x, sig_SNPs, race, params), 
                       axis=1, result_type='expand')

    # Reshape data into the right structure
    geno = geno.stack().reset_index().rename({"level_0": "gene", 
                                              "level_1": "nwd_id", 0: "genotype"}, axis=1)
    pheno = pheno.stack().reset_index().rename({"level_0": "gene", 
                                                "level_1": "nwd_id", 0: "expression"}, axis=1)
    merged_simulated_data = pd.merge(geno, pheno)
    merged_simulated_data = pd.merge(merged_simulated_data, sig_SNPs, 
                                     left_on="gene", right_index=True)
    merged_simulated_data = pd.merge(merged_simulated_data, race, 
                                     left_on="nwd_id", right_index=True)
    merged_simulated_data[["genotype_Eur", "genotype_Afr"]] = \
        merged_simulated_data.apply(lambda x: pd.Series([int(x.genotype[0]), int(x.genotype[2])]), axis=1)
    merged_simulated_data["intercept_Eur"] = None
    merged_simulated_data["intercept_Afr"] = None

    # Write simulated data to file
    merged_simulated_data.to_csv(args.out, sep='\t', index=False)

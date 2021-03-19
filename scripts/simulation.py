import pandas as pd
import numpy as np
import argparse
import scipy.stats
from identify_hits import FDR_threshold

# Specify parameters for the beta distribution used to sample global ancestry
SHAPE_1 = 7.9
SHAPE_2 = 2.1

def create_MAF_dist(afr_hit_path, eur_hit_path):
    # Create MAF distribution from empirical dataset
    afr_hits = pd.read_csv(afr_hit_path, sep='\t')
    eur_hits = pd.read_csv(eur_hit_path, sep='\t')
    hits = pd.merge(afr_hits, eur_hits, on=["Gene", "ID"], 
                    suffixes = ("_Afr", "_Eur"))
    hits = np.array(hits[["maf_Afr", "maf_Eur"]])
    return(hits)

def simulate_tag_SNP(n, maf_dist, params):
    """Simulate a mix of causal/trait-associated/significant SNPs and 
       non-causal/non-associated/insignificant SNPs."""
    # Simulate betas for trait-associated SNPs as multivariate normal
    beta_cov_matrix = np.array([[params["var_afr_beta"], params["covariance"]], 
                                [params["covariance"], params["var_eur_beta"]]])
    beta_mean = np.array([params["mean_afr_beta"], params["mean_eur_beta"]])
    tag_SNPs = pd.DataFrame(np.random.multivariate_normal(beta_mean, beta_cov_matrix, n), 
                            columns=["Beta_Afr", "Beta_Eur"])

    # Simulate non-associated SNPs
    tag_SNPs.loc[0:int(n * params["frac_non_causal"]), ["Beta_Afr", "Beta_Eur"]] = 0

    # Simulate MAF by uniformly randomly sampling from empirical distribution
    # TODO: we sample MAF for associated AND non-associated SNPs from the
    # empirical distribution for associated SNPs - is this a problem?
    tag_SNPs["MAF_Afr"], tag_SNPs["MAF_Eur"] = None, None
    tag_SNPs[["MAF_Afr", "MAF_Eur"]] = maf_dist[np.random.randint(maf_dist.shape[0], size=n), :]

    return(tag_SNPs)

def simulate_estimated_beta(row, n_idv, error, pop):
    """Simulate the estimated beta for a tag SNP from a normal distribution,
       with mean equal to the true beta and variance determined by the
       pre-specified error parameter."""
    var_X = 2 * row["MAF_" + pop] * (1 - row["MAF_" + pop])
    sd_beta = np.sqrt(error / (n_idv * var_X))
    est_beta = np.random.normal(row["Beta_" + pop], sd_beta)
    return(est_beta) 

def simulate_beta_sample_var(row, n_idv, error, pop):
    """Simulate residuals from normal distribution based on pre-specified
       error parameter. Simulate genotypes from binomial distribution based on
       MAF of SNP. Use simulated residuals and genotypes to simulate sample 
       variance of beta."""
    sim_resid = np.random.normal(0, np.sqrt(error), size=n_idv)
    sim_geno = np.random.binomial(2, row["MAF_" + pop], size=n_idv)
    sample_var = (np.sum(sim_resid ** 2)) / ((n_idv - 2) * (n_idv * np.var(sim_geno)))
    return(sample_var)

def simulate_estimation(tag_SNPs, group_name, n_idv, error):
    """Simulate estimated betas and sample variance of all SNP betas. Compute 
       associated T-statistics and p-values."""
    if group_name == "Afr":
        pop = group_name
    else:
        pop = "Eur"
    tag_SNPs["effect_" + group_name] = tag_SNPs.apply(lambda row: 
        simulate_estimated_beta(row, n_idv, error, pop), axis=1)
    tag_SNPs["Sample_Var_" + group_name] = tag_SNPs.apply(lambda row: 
        simulate_beta_sample_var(row, n_idv, error, pop), axis=1)
    tag_SNPs["Tval_" + group_name] = tag_SNPs.apply(lambda row: 
        row["effect_" + group_name] / np.sqrt(row["Sample_Var_" + group_name]), axis=1)
    tag_SNPs["Pval_" + group_name] = tag_SNPs.apply(lambda row: 
        1 - scipy.stats.t.cdf(row["Tval_" + group_name], df=2 * n_idv - 2) if 
        row["Tval_" + group_name] > 0 else 
        scipy.stats.t.cdf(row["Tval_" + group_name], df=2 * n_idv - 2), axis=1)
    return(tag_SNPs)

def simulate_significant_SNPs(maf_dist, params):
    """Simulate set of significant SNPs by first simulating a mixture of 
       associated SNPs and non-associated SNPs, simulating estimation procedure,
       and repeating until correct number of associated/significant SNPs are
       obtained."""
    tag_SNPs = simulate_tag_SNP(1500, maf_dist, params)
    tag_SNPs = simulate_estimation(tag_SNPs, "Asc", 
                                   params["n_asc_idv"], params["var_error_eur"])
    sig_SNPs = FDR_threshold(tag_SNPs, "Pval_Asc", params["fdr_threshold"])

    # Continue simulating tag SNPs until we obtain the correct number of significant SNPs
    while sig_SNPs.shape[0] < params["n_sig_snp"]:
        new_SNPs = simulate_tag_SNP(100, maf_dist, params)
        new_SNPs = simulate_estimation(new_SNPs, "Asc", 
                                       params["n_asc_idv"], params["var_error_eur"])
        tag_SNPs = pd.concat([tag_SNPs, new_SNPs])
        sig_SNPs = FDR_threshold(tag_SNPs, "Pval_Asc", params["fdr_threshold"])
    return(sig_SNPs)

def simulate_ganc(row):
    """Simulate global ancestry from beta distribution specified by shape
       parameters at top of script."""
    ganc = row[0] * np.random.beta(SHAPE_1, SHAPE_2)
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
        ind_race = race_df.loc[ind, "Race_AA"]
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
    parser.add_argument('--var_afr', help="specified variance of African tag SNP effects", default=None)
    parser.add_argument('--var_eur', help="specified variance of European tag SNP effects", default=None)
    parser.add_argument('--correlation', help="specified correlation between African and European tag SNP effects", default=None)
    parser.add_argument('--n_snp', help="specified number of significant SNPs to simulate", default=None)
    parser.add_argument('--n_idv', help="comma-delimited number of ascertainment and reestimation individuals to simulate", default=None)
    parser.add_argument('--error_afr', help="specified error term", default=None)
    parser.add_argument('--error_eur', help="specified error term", default=None)
    parser.add_argument('--delta', help="specified delta term", default=None)
    parser.add_argument('--afr_hits', help="African hits")
    parser.add_argument('--eur_hits', help="European hits")
    parser.add_argument('--out', help="output file")
    args = parser.parse_args()

    # Parse command-line arguments and specify parameters for simulation
    params = {"mean_afr_beta": 0, "mean_eur_beta": 0, "fdr_threshold": .05, "frac_non_causal": 0.8}
    if args.var_afr:
        params["var_afr_beta"] = float(args.var_afr)
    else:
        params["var_afr_beta"] = 0.17
    if args.var_eur:
        params["var_eur_beta"] = float(args.var_eur)
    else:
        params["var_eur_beta"] = 0.25
    if args.correlation:
        params["covariance"] = float(args.correlation) * ((params["var_afr_beta"] * params["var_eur_beta"]) ** 0.5)
    else:
        params["covariance"] = 0.17
    if args.n_snp:
        params["n_sig_snp"] = int(args.n_snp)
    else:
        params["n_sig_snp"] = 3200
    if args.n_idv:
        n_idv = args.n_idv.strip().split('_')
        params["n_asc_idv"] = int(n_idv[0])
        params["n_est_idv"] = int(n_idv[1])
    else:
        params["n_asc_idv"] = 200
        params["n_est_idv"] = 150
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

    # Create MAF distribution from input data
    maf_dist = create_MAF_dist(args.afr_hits, args.eur_hits)

    # Simulate African and European betas for significant SNPs, with 
    # ascertainment performed in Europeans
    sig_SNPs = simulate_significant_SNPs(maf_dist, params)
    sig_SNPs = simulate_estimation(sig_SNPs, "Eur", params["n_est_idv"], params["var_error_eur"])
    sig_SNPs = simulate_estimation(sig_SNPs, "Afr", params["n_est_idv"], params["var_error_afr"])
    n_sig_SNPs = sig_SNPs.shape[0]

    # Simulate local ancestry, genotypes, and phenotypes for each of the 
    # significant SNPs
    race = pd.DataFrame([anc for anc in [0, 1] for _ in range(params["n_est_idv"])], 
                        columns=["Race_AA"])
    ganc = pd.DataFrame(race.apply(simulate_ganc, axis=1))
    lanc = ganc.apply(lambda x: simulate_lanc(x, n_sig_SNPs), axis=1, 
                      result_type='expand').T
    geno = lanc.apply(lambda x: simulate_genotypes(x, sig_SNPs), axis=1, 
                      result_type='expand')
    pheno = geno.apply(lambda x: simulate_phenotypes(x, sig_SNPs, race, params), 
                       axis=1, result_type='expand')

    # Reshape data into the right structure
    geno = geno.stack().reset_index().rename({"level_0": "Gene", 
                                              "level_1": "Ind", 0: "Genotype"}, axis=1)
    pheno = pheno.stack().reset_index().rename({"level_0": "Gene", 
                                                "level_1": "Ind", 0: "Expression"}, axis=1)
    merged_simulated_data = pd.merge(geno, pheno)
    merged_simulated_data = pd.merge(merged_simulated_data, sig_SNPs, 
                                     left_on="Gene", right_index=True)
    merged_simulated_data = pd.merge(merged_simulated_data, race, 
                                     left_on="Ind", right_index=True)
    merged_simulated_data[["Genotype_Eur", "Genotype_Afr"]] = \
        merged_simulated_data.apply(lambda x: pd.Series([int(x.Genotype[0]), int(x.Genotype[2])]), axis=1)
    merged_simulated_data["intercept_Eur"] = None
    merged_simulated_data["intercept_Afr"] = None

    # Write simulated data to file
    merged_simulated_data.to_csv(args.out, sep='\t', index=False)

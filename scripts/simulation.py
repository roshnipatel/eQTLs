import pandas as pd
import numpy as np
import argparse
import scipy.stats
from identify_hits import FDR_threshold

# Specify constants
VAR_AFR_BETA = .3
VAR_EUR_BETA = .3
MEAN_AFR_BETA = 0
MEAN_EUR_BETA = 0
N_ASC_IDV = 200
N_EST_IDV = 150
N_SIG_SNP = 1500
FDR_THRESHOLD = .1
FRAC_NON_CAUSAL = .8

# Specify parameters for the beta distribution used to sample global ancestry
SHAPE_1 = 7.9
SHAPE_2 = 2.1

def create_MAF_dist(afr_hit_path, eur_hit_path):
    # Create MAF distribution from empirical dataset
    afr_hits = pd.read_csv(afr_hit_path, sep='\t')
    eur_hits = pd.read_csv(eur_hit_path, sep='\t')
    hits = pd.merge(afr_hits, eur_hits, on=["Gene", "ID"], suffixes = ("_Afr", "_Eur"))
    hits = np.array(hits[["maf_Afr", "maf_Eur"]])
    return(hits)

def simulate_tag_SNP(n, maf_dist):
    # Simulate betas (effect sizes) as multivariate normal
    beta_cov_matrix = np.array([[VAR_AFR_BETA, COVARIANCE], [COVARIANCE, VAR_EUR_BETA]])
    beta_mean = np.array([MEAN_AFR_BETA, MEAN_EUR_BETA])
    tag_SNPs = pd.DataFrame(np.random.multivariate_normal(beta_mean, beta_cov_matrix, n), columns=["Beta_Afr", "Beta_Eur"])
    # Simulate non-causal variants
    tag_SNPs.loc[0:int(n * FRAC_NON_CAUSAL), ["Beta_Afr", "Beta_Eur"]] = 0
    # Simulate MAF by uniformly randomly sampling from empirical distribution
    tag_SNPs["MAF_Afr"], tag_SNPs["MAF_Eur"] = None, None
    tag_SNPs[["MAF_Afr", "MAF_Eur"]] = maf_dist[np.random.randint(maf_dist.shape[0], size=n), :]
    return(tag_SNPs)

def simulate_estimated_beta(row, n_idv, error, pop):
    # Simulate estimated beta as normal (true variance is "known") 
    var_X = 2 * row["MAF_" + pop] * (1 - row["MAF_" + pop])
    sd_beta = np.sqrt(error / (n_idv * var_X))
    est_beta = np.random.normal(row["Beta_" + pop], sd_beta)
    return(est_beta) 

def simulate_beta_sample_var(row, n_idv, error, pop):
    # Simulate residuals and genotypes (based on MAF of variant) to simulate sample variance of beta
    sim_resid = np.random.normal(0, np.sqrt(error), size=n_idv)
    sim_geno = np.random.binomial(2, row["MAF_" + pop], size=n_idv)
    sample_var = (np.sum(sim_resid ** 2)) / ((n_idv - 2) * (n_idv * np.var(sim_geno)))
    return(sample_var)

def simulate_estimation(tag_SNPs, group_name, n_idv, error):
    # Simulate estimated betas and sample variance of betas. Compute associated T-statistics and p-values.
    if group_name == "Afr":
        pop = group_name
    else:
        pop = "Eur"
    tag_SNPs["effect_" + group_name] = tag_SNPs.apply(lambda row: simulate_estimated_beta(row, n_idv, error, pop), axis=1)
    tag_SNPs["Sample_Var_" + group_name] = tag_SNPs.apply(lambda row: simulate_beta_sample_var(row, n_idv, error, pop), axis=1)
    tag_SNPs["Tval_" + group_name] = tag_SNPs.apply(lambda row: row["effect_" + group_name] / np.sqrt(row["Sample_Var_" + group_name]), axis=1)
    tag_SNPs["Pval_" + group_name] = tag_SNPs.apply(lambda row: 1 - scipy.stats.t.cdf(row["Tval_" + group_name], df=2 * n_idv - 2) if row["Tval_" + group_name] > 0 else scipy.stats.t.cdf(row["Tval_" + group_name], df=2 * n_idv - 2), axis=1)
    return(tag_SNPs)

def FDR_threshold(tag_SNPs, colname, fdr):
    n_SNPs = tag_SNPs.shape[0]
    def check_significance(row):
        pval_threshold = (row.name + 1) / n_SNPs * fdr
        return(not (row[colname] > pval_threshold))
    ordered_SNPs = tag_SNPs.sort_values(by=colname, ignore_index=True)
    ordered_SNPs["Significant"] = ordered_SNPs.apply(check_significance, axis=1)
    sig_SNPs = ordered_SNPs[ordered_SNPs.Significant == True]
    return(sig_SNPs)

def simulate_significant_SNPs(maf_dist):
    global total_sim_SNPs
    total_sim_SNPs = 1500
    tag_SNPs = simulate_tag_SNP(1500, maf_dist)
    tag_SNPs = simulate_estimation(tag_SNPs, "Asc", N_ASC_IDV, VAR_ERROR)
    sig_SNPs = FDR_threshold(tag_SNPs, "Pval_Asc", FDR_THRESHOLD)
    # Continue simulating tag SNPs until we obtain the correct number of significant SNPs
    while sig_SNPs.shape[0] < N_SIG_SNP:
        new_SNPs = simulate_tag_SNP(100, maf_dist)
        new_SNPs = simulate_estimation(new_SNPs, "Asc", N_ASC_IDV, VAR_ERROR)
        tag_SNPs = pd.concat([tag_SNPs, new_SNPs])
        sig_SNPs = FDR_threshold(tag_SNPs, "Pval_Asc", FDR_THRESHOLD)
        total_sim_SNPs += 100
    return(sig_SNPs)

def simulate_ganc(row):
    ganc = row[0] * np.random.beta(SHAPE_1, SHAPE_2)
    return(ganc)

def simulate_lanc(row, n):
    # Simulate local ancestry as binomial conditional on race/population
    lanc = np.random.binomial(2, row[0], size=n)
    return(lanc)

def simulate_genotypes(row, SNP_df):
    geno = []
    # Simulate genotypes as binomial conditional on local ancestry
    for lanc in row:
        geno_Afr = np.random.binomial(lanc, SNP_df.loc[row.name, "MAF_Afr"])
        geno_Eur = np.random.binomial(2 - lanc, SNP_df.loc[row.name, "MAF_Eur"])
        geno.append(str(geno_Afr) + "-" + str(geno_Eur))
    return(geno)

def simulate_phenotypes(row, SNP_df):
    pheno = []
    # Simulate phenotypes as betas weighted by genotypes, plus random, normally-distributed noise
    for geno in row:
        pheno_Afr = int(geno[0]) * SNP_df.loc[row.name, "Beta_Afr"]
        pheno_Eur = int(geno[2]) * SNP_df.loc[row.name, "Beta_Eur"]
        noise = np.random.normal(0, np.sqrt(VAR_ERROR))
        pheno.append(pheno_Afr + pheno_Eur + noise)
    return(pheno)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--covariance', help="specified covariance between African and European tag SNP effects")
    parser.add_argument('--error', help="specified error term")
    parser.add_argument('--afr_hits', help="African hits")
    parser.add_argument('--eur_hits', help="European hits")
    parser.add_argument('--out', help="output file")
    args = parser.parse_args()

    # Parse constants
    global COVARIANCE
    COVARIANCE = float(args.covariance)
    global VAR_ERROR
    VAR_ERROR = float(args.error)

    # Create MAF distribution from input data
    maf_dist = create_MAF_dist(args.afr_hits, args.eur_hits)

    # Simulate African and European betas for significant SNPs, with ascertainment performed in Europeans
    sig_SNPs = simulate_significant_SNPs(maf_dist)
    sig_SNPs = simulate_estimation(sig_SNPs, "Eur", N_EST_IDV, VAR_ERROR)
    sig_SNPs = simulate_estimation(sig_SNPs, "Afr", N_EST_IDV, VAR_ERROR)
    n_sig_SNPs = sig_SNPs.shape[0]

    # Simulate local ancestry, genotypes, and phenotypes for each of the significant SNPs/genes
    race = pd.DataFrame([anc for anc in [0, 1] for _ in range(N_EST_IDV)], columns=["Race_AA"])
    ganc = pd.DataFrame(race.apply(simulate_ganc, axis=1))
    lanc = ganc.apply(lambda x: simulate_lanc(x, n_sig_SNPs), axis=1, result_type='expand').T
    geno = lanc.apply(lambda x: simulate_genotypes(x, sig_SNPs), axis=1, result_type='expand')
    pheno = geno.apply(lambda x: simulate_phenotypes(x, sig_SNPs), axis=1, result_type='expand')

    # Reshape data into the right structure
    geno = geno.stack().reset_index().rename({"level_0": "Gene", "level_1": "Ind", 0: "Genotype"}, axis=1)
    pheno = pheno.stack().reset_index().rename({"level_0": "Gene", "level_1": "Ind", 0: "Expression"}, axis=1)
    merged_simulated_data = pd.merge(geno, pheno)
    merged_simulated_data = pd.merge(merged_simulated_data, sig_SNPs, left_on="Gene", right_index=True)
    merged_simulated_data = pd.merge(merged_simulated_data, race, left_on="Ind", right_index=True)
    merged_simulated_data[["Genotype_Afr", "Genotype_Eur"]] = merged_simulated_data.apply(lambda x: pd.Series([int(x.Genotype[0]), int(x.Genotype[2])]), axis=1)
    merged_simulated_data["intercept_Eur"] = None
    merged_simulated_data["intercept_Afr"] = None

    # Write simulated data to file
    merged_simulated_data.to_csv(args.out, sep='\t', index=False)

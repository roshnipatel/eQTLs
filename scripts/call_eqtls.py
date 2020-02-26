import sys
import pandas as pd
import argparse
import numpy as np
import statsmodels.api as sm

# TODO add MAF threshold; include MAF in final output
# TODO handle covariates

regression_stat_cols = ["effect", "se", "r2", "pval"]
geno_info_cols = ["ID", "ma_samp", "maf"]

def get_ind(sample_filepath):
    # Read sample IDs from file into list
    sample_file = open(sample_filepath, 'r')
    ind = [i.strip() for i in sample_file.readlines()]
    sample_file.close()
    return(ind)

def parse_expression(ind):
    # Parse expression data input and filter for desired sample IDs
    for line in sys.stdin:
        if line[0] == '#':
            header = line.strip().split()
        else:
            exp = line.strip().split()
    exp = pd.DataFrame(exp, columns=header)
    exp = exp[ind]
    return(exp)

def parse_covariates(cov_file_path, ind):
    cov = pd.read_csv(cov_file_path, delimiter='\t')
    cov = cov[ind]
    return(cov)

def parse_genotypes(vcf_path, ind):
    # Parse genotype data input and filter for desired sample IDs
    geno_list = []
    with open(vcf_path, 'r') as geno_file:
        for line in geno_file:
            if line[:2] == "##":
                continue
            if line[0] == "#":
                header = line.split()
            else:
                vcf_geno = line.split()
                int_geno = [int(gt[0]) + int(gt[2]) for gt in vcf_geno] # Convert genotype strings in VCF file into integer allele counts
                geno_list.append(int_geno)
    if not geno_list: # Checks for the case in which there are no SNPs within 100 kb of the given gene's TSS
        return(None)
    else:
        geno = pd.DataFrame(geno_list, columns=header)["ID" + ind] # Drop all columns except for SNP ID and individuals in sample file
        return(geno)

def maf_filter(geno, ind, maf_thresh, ma_samp_thresh):
    alleles = geno[ind] # Create separate DataFrame to store allele counts
    
    # Record number of individuals with minor allele for each SNP
    geno["ma_samp"] = alleles[alleles != 0].count(axis=1)
    
    # Record minor allele frequency for each SNP
    n_ind = len(ind)
    geno["maf"] = alleles.sum(axis=1).divide(n_ind)

    # Remove rows that do not meet both minor allele thresholds
    geno = geno[(geno.ma_samp >= ma_samp_thresh) & (geno.maf >= maf_thresh)]    
    return(geno)

def perform_regression(geno, exp, ind):
    # Regress genotypes onto expression data
    def marginal_test(row):
        model = sm.OLS(exp.iloc[0], row[ind])
        results = model.fit()
        effect = results.params[0]
        se = results.bse[0]
        r2 = results.rsquared
        pval = results.pvalues[0]
        return(effect, se, pval, r2)

    res = pd.DataFrame()
    res[regression_stat_cols] = geno.apply(marginal_test, axis=1, result_type='expand')
    res[geno_info_cols] = geno[geno_info_cols]
    return(res)

def perm_adjust_pval(n, min_pval, geno, exp, ind):
    def permutation(geno, exp, ind):
        randomized_ind = np.random.permutation(ind)
        exp = exp[randomized_ind]
        perm_min_pval = perform_regression(geno, exp, ind)["pval"].min(axis=1)
        return(perm_min_pval)

    null_pval_dist = [permutation(geno, exp, ind) for _ in range(n)]
    adj_pval = (null_pval_dist[null_pval_dist < min_pval].count().sum() + 1) / (n + 1)
    return(adj_pval)

def correct_covariates(exp, ind, cov):
    model = sm.OLS(exp.iloc[0], cov)
    results = model.fit()
    resid_exp = results.resid
    return(resid_exp)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genotypes', help='Genotypes of all individuals for SNPs within 100 kb of gene TSS.')
    parser.add_argument('--samples', help='File with sample IDs for individuals to use when calling eQTLs.')
    parser.add_argument('--n_perm', help='Number of permutations to run.')
    parser.add_argument('--covariates', default=None, help='Covariates to include in regression.')
    parser.add_argument('--maf', default=0.05, help='MAF threshold for estimating variant effect size.')
    parser.add_argument('--ma_samples', default=5, help='Minor allele sample count threshold for estimating variant effect size.')
    parser.add_argument('--out', help='Output file name.')
    args = parser.parse_args()
    
    # Record number of permutations to run
    try:
        n_perm = int(args.n_perm)
    except:
        print("Must specify integer number of permutations to run.")

    ind_IDs = get_ind(args.samples)

    expression = parse_expression(ind_IDs)
    if args.covariates is not None:
        covariates = parse_covariates(args.covariates, ind_IDs)
        expression = correct_covariates(expression, ind_IDs, covariates)

    genotypes = parse_genotypes(args.genotypes)
    if genotypes is not None:
        genotypes = maf_filter(genotypes, ind_IDs, args.maf, args.ma_samples)
        results = perform_regression(genotypes, expression, ind_IDs)
        candidate = results[results.pval == results.pval.min()]
    
        cand_pval = candidate["pval"].iloc[0]
        candidate["perm_pval"] = perm_adjust_pval(n_perm, cand_pval, genotypes, expression, ind_IDs)
        candidate.to_csv(args.out, sep='\t')

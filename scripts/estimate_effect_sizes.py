import sys
from datetime import datetime
import os
import pandas as pd
import argparse
import gzip
import numpy as np
import statsmodels.api as sm

def get_ind(sample_filepath):
    """Read sample IDs from file into list."""
    sample_file = open(sample_filepath, 'r')
    ind = [i.strip() for i in sample_file.readlines()]
    sample_file.close()
    return(ind)

def parse_expression(exp_path, ind):
    """Parse expression data input and filter for desired sample IDs."""
    exp = pd.read_csv(exp_path, sep='\t')
    exp = exp[ind]
    exp = exp.dropna(axis=1)
    filtered_ind = exp.columns
    exp = exp.T.squeeze().astype(float)
    return(exp, filtered_ind)

def parse_covariates(cov_file_path, ind):
    """Parse covariates and filter for desired sample IDs."""
    cov = pd.read_csv(cov_file_path, delimiter='\t', index_col="nwd_id")
    cov = cov.loc[ind]
    return(cov)

def parse_genotypes(vcf_path, ind=None, concat=True):
    """Parse genotype data input and filter for desired sample IDs."""
    geno_list = []
    with gzip.open(vcf_path, 'rt') as geno_file:
        for line in geno_file:
            if line[:2] == "##":
                continue
            if line[0] == "#":
                line = line.split()
                header = [line[2]] + line[9:]
            else:
                line = line.split()
                ID = line[2]
                vcf_geno = line[9:]
                if concat:
                    # Convert genotype strings in VCF file into integer counts
                    int_geno = [int(gt[0]) + int(gt[2]) for gt in vcf_geno] 
                    geno_list.append([ID] + int_geno)
                else:
                    geno_list.append([ID] + vcf_geno)
    geno = pd.DataFrame(geno_list, columns=header).set_index(header[0])
    if ind is not None:
        geno = geno[ind]
    return(geno)

def maf_filter(geno, maf_thresh, ma_samp_thresh):
    """Exclude genotypes that do not meet minimum maf and minor allele sample
       thresholds."""
    n_ind = geno.columns.shape[0]

    # Record alternate allele count and frequency
    maf_table = pd.DataFrame({"alt_samp": geno[geno != 0].count(axis=1),
                              "alt_af": geno.sum(axis=1).divide(2 * n_ind)})

    # Record minor allele frequency and count
    maf_table["ma_samp"] = maf_table.apply(lambda row: min(row.alt_samp, n_ind - row.alt_samp), axis=1)
    maf_table["maf"] = maf_table.apply(lambda row: min(row.alt_af, 1 - row.alt_af), axis=1)
     
    # Remove SNPs that do not meet both minor allele thresholds
    maf_table = maf_table[(maf_table.ma_samp >= ma_samp_thresh) & 
                          (maf_table.maf >= maf_thresh)]

    # Apply transformation to table of genotypes
    geno = geno.loc[maf_table.index]

    return(geno, maf_table)

def regress_genotypes(geno, exp):
    """Regress genotypes onto expression data."""
    def marginal_test(x, y):
        x = pd.DataFrame(x.astype(float))
        x["intercept"] = 1
        model = sm.OLS(y, x)
        results = model.fit()
        effect = results.params[0]
        se = results.bse[0]
        r2 = results.rsquared
        pval = results.pvalues[0]
        return(effect, se, r2, pval)

    regression_stat_cols = ["effect", "se", "r2", "pval"]
    res = pd.DataFrame(columns=regression_stat_cols, index=geno.index)
    res[regression_stat_cols] = geno.apply(lambda row: 
        marginal_test(row, exp), axis=1, result_type='expand')
    return(res)

def perm_adjust_pval(n, min_pval, geno, exp, ind):
    """Adjust p-values using permutations of genotype/expression data. Not
       currently used in pipeline."""
    def permutation(geno, exp, ind):
        # Randomize data, decoupling individuals' expression from genotypes
        randomized_ind = np.random.permutation(ind)
        exp = exp[randomized_ind]
        geno.columns = randomized_ind

        # Regress randomized expression on genotypes
        perm_min_pval = regress_genotypes(geno, exp)["pval"].min()
        return(perm_min_pval)

    null_pval_dist = pd.Series([permutation(geno, exp, ind) for _ in range(n)])
    adj_pval = (null_pval_dist[null_pval_dist < min_pval].count().sum() + 1) / n
    return(adj_pval)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genotypes')
    parser.add_argument('--phenotypes')
    parser.add_argument('--samples', help='File with sample IDs')
    parser.add_argument('--covariates', default=None)
    parser.add_argument('--maf', default=0.05)
    parser.add_argument('--ma_samples', default=5)
    parser.add_argument('--n_perm', default=0, help='Number of permutations to run.')
    parser.add_argument('--out')
    args = parser.parse_args()

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("Start Time =", current_time)

    # If regression fails, set command to create blank output file
    fail_command = "touch {0}".format(args.out)

    # Record number of permutations to run
    n_perm = int(args.n_perm)
    
    # Load expression and genotype data for ALL individuals, then select the
    # individuals for the current regression (ind_IDs)
    ind_IDs = get_ind(args.samples)
    expression, ind_IDs = parse_expression(args.phenotypes, ind_IDs)
    genotypes = parse_genotypes(args.genotypes, ind_IDs)
    genotypes, maf_table = maf_filter(genotypes, args.maf, args.ma_samples)

    # Exit if there are no SNPs within 100 kb of the gene's TSS,
    # creating blank output file so that Snakemake doesn't throw a fit
    if genotypes.empty: 
        os.system(fail_command)
        sys.exit()

    # Residualize expression on covariates
    if args.covariates is not None:
        covariates = parse_covariates(args.covariates, ind_IDs)
        covariates["intercept"] = 1
        model = sm.OLS(expression, covariates)
        results = model.fit()
        expression = results.resid

    # Regress expression data on genotype data 
    results = regress_genotypes(genotypes, expression)
    results = pd.merge(results, maf_table, left_index=True, right_index=True)

    # Adjust p-values with permutations if desired and write results to output
    if n_perm == 0:
        results.to_csv(args.out, sep='\t')
    else:
        # Select SNP with min pval (randomly choose in case of tie)
        candidate = results[results.pval == results.pval.min()].sample(1).T.squeeze()
        cand_pval = candidate.pval

        # If candidate SNP is significant, perform permutations
        if cand_pval < .05:
            perm_pval = perm_adjust_pval(n_perm, cand_pval, 
                                         genotypes, expression, ind_IDs)
        else:
            perm_pval = 1
            n_perm = 0
        candidate["perm_pval"] = perm_pval
        candidate["n_perm"] = n_perm
        candidate.to_frame().T.to_csv(args.out, sep='\t')

    now = datetime.now()
    current_time = now.strftime("%H:%M:%S")
    print("End Time =", current_time)

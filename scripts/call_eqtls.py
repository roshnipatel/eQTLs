import sys
import os
import pandas as pd
import argparse
import gzip
import numpy as np
import statsmodels.api as sm

regression_stat_cols = ["effect", "se", "r2", "pval"]
geno_info_cols = ["ma_samp", "maf"]

def get_ind(sample_filepath):
    # Read sample IDs from file into list
    sample_file = open(sample_filepath, 'r')
    ind = [i.strip() for i in sample_file.readlines()]
    sample_file.close()
    return(ind)

def parse_expression(exp_path, ind):
    # Parse expression data input and filter for desired sample IDs
    exp = pd.read_csv(exp_path, sep='\t')
    exp = exp[ind]
    # exp = exp - exp.mean(axis=1).iloc[0] # Mean-center expression data
    return(exp)

def parse_covariates(cov_file_path, ind):
    cov = pd.read_csv(cov_file_path, delimiter='\t', index_col="NWDID").T
    cov = cov[ind]
    return(cov)

def parse_genotypes(vcf_path, ind=None, concat=True):
    # Parse genotype data input and filter for desired sample IDs
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
                    int_geno = [int(gt[0]) + int(gt[2]) for gt in vcf_geno] # Convert genotype strings in VCF file into integer allele counts
                    geno_list.append([ID] + int_geno)
                else:
                    geno_list.append([ID] + vcf_geno)
    geno = pd.DataFrame(geno_list, columns=header).set_index(header[0])
    if ind is not None:
        geno = geno[ind]
    return(geno)

def maf_filter(geno, ind, maf_thresh, ma_samp_thresh):
    alleles = geno[ind] # Create separate DataFrame to store allele counts
    
    # Record number of individuals with alternate allele for each SNP
    geno["alt_samp"] = alleles[alleles != 0].count(axis=1)
    
    # Record alternate allele frequency for each SNP
    n_ind = len(ind)
    geno["alt_af"] = alleles.sum(axis=1).divide(2 * n_ind)

    # Record minor allele frequency and count
    geno["ma_samp"] = pd.DataFrame([geno["alt_samp"], n_ind - geno["alt_samp"]]).T.min(axis=1)
    geno["maf"] = pd.DataFrame([geno["alt_af"], 1 - geno["alt_af"]]).T.min(axis=1)

    # Remove rows that do not meet both minor allele thresholds
    geno = geno[(geno.ma_samp >= ma_samp_thresh) & (geno.maf >= maf_thresh)][list(ind) + geno_info_cols]
    return(geno)

def perform_regression(geno, exp, ind, cov=None):
    # Regress genotypes onto expression data
    def marginal_test(row, exp, cov=None):
        y = exp.iloc[0].astype(float)
        x = pd.DataFrame(row.astype(float))
        x["int"] = 1
        if cov is not None:
            x = pd.merge(x, cov.T, left_index=True, right_index=True)
        model = sm.OLS(y, x)
        results = model.fit()
        effect = results.params[0]
        se = results.bse[0]
        r2 = results.rsquared
        pval = results.pvalues[0]
        return(effect, se, r2, pval)
    res = pd.DataFrame(columns=regression_stat_cols, index=geno.index)
    res[regression_stat_cols] = geno.apply(lambda x: marginal_test(x, exp, cov), axis=1, result_type='expand')
    return(res)

def perm_adjust_pval(n, min_pval, geno, exp, ind):
    def permutation(geno, exp, ind):
        randomized_ind = np.random.permutation(ind)
        n_ind = len(ind)
        dummy_ind = list(range(n_ind))
        exp = exp[randomized_ind]
        geno = geno[ind]
        exp.columns, geno.columns = dummy_ind, dummy_ind # Need to match column names in order to call sm.OLS
        perm_min_pval = perform_regression(geno, exp, dummy_ind)["pval"].min()
        return(perm_min_pval)
    null_pval_dist = pd.Series([permutation(geno, exp, ind) for _ in range(n)])
    adj_pval = (null_pval_dist[null_pval_dist < min_pval].count().sum() + 1) / (n + 1)
    return(adj_pval)

def remove_nan(geno, exp):
    all_ind = exp.columns
    exp_dropna = exp.dropna(axis=1)
    ind_dropna = exp_dropna.columns
    dropped_ind = all_ind.difference(ind_dropna)
    geno_dropna = geno.drop(dropped_ind, axis=1)
    return(geno_dropna, exp_dropna, ind_dropna)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--genotypes', help='Genotypes of all individuals for SNPs within 100 kb of gene TSS.')
    parser.add_argument('--phenotypes', help='Phenotypes of all individuals for the given gene.')
    parser.add_argument('--samples', help='File with sample IDs for individuals to use when calling eQTLs.')
    parser.add_argument('--n_perm', default=1000, help='Number of permutations to run.')
    parser.add_argument('--type', help='est (estimation; no permutations) or asc (ascertainment; permutations performed)')
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
    
    fail_command = "touch {0}".format(args.out)
    
    print("Reading input files.")
    ind_IDs = get_ind(args.samples)
    expression = parse_expression(args.phenotypes, ind_IDs)
    genotypes = parse_genotypes(args.genotypes, ind_IDs)
    
    print("Done reading input files.")
    if genotypes is None: # Checks for the case in which there are no SNPs within 100 kb of the given gene's TSS
        os.system(fail_command)
    else: 
        genotypes, expression, ind_IDs = remove_nan(genotypes, expression)
        genotypes = maf_filter(genotypes, ind_IDs, args.maf, args.ma_samples)
        if genotypes.empty: # Checks for the case in which there are no SNPs in the desired MAF window
            os.system(fail_command)
        else:
            if args.covariates is not None:
                covariates = parse_covariates(args.covariates, ind_IDs)
                results = perform_regression(genotypes, expression, ind_IDs, covariates)
            else:
                results = perform_regression(genotypes, expression, ind_IDs)
            results[geno_info_cols] = genotypes[geno_info_cols]
            results.to_csv(args.out, sep='\t')

#         print("Finished performing regression.")
#         if args.type == 'ascertainment':
#             print("Estimating eQTL for ascertainment dataset.")
#             candidate = results[results.pval == results.pval.min()].sample(1) # Randomly choose SNP with min pval in case of tie
#             cand_pval = candidate["pval"].iloc[0]
#             if cand_pval < .10:
#                 print("Adjusting p-value; performing 100 permutations.")
#                 perm_pval = perm_adjust_pval(n_perm, cand_pval, genotypes, expression, ind_IDs)
#             else:
#                 print("No permutations performed; initial p-value > .10.")
#                 candidate[["perm_pval", "n_perm"]] = [None, 0]
#                 candidate.to_csv(args.out, sep='\t')
#             # while perm_pval < .10 and n_perm <= 1000:
#             #     n_perm = n_perm * 10
#             #     print("Adjusting p-value; performing 1000 permutations.")
#             #     perm_pval = perm_adjust_pval(n_perm, cand_pval, genotypes, expression, ind_IDs)
#             candidate["perm_pval"] = perm_pval
#             candidate["n_perm"] = n_perm
#             candidate.to_csv(args.out, sep='\t')
#         else:
#             results.to_csv(args.out, sep='\t')

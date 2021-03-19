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
    """Read sample IDs from file into list."""
    sample_file = open(sample_filepath, 'r')
    ind = [i.strip() for i in sample_file.readlines()]
    sample_file.close()
    return(ind)

def parse_expression(exp_path, ind):
    """Parse expression data input and filter for desired sample IDs."""
    exp = pd.read_csv(exp_path, sep='\t')
    exp = exp[ind]
    return(exp)

def parse_covariates(cov_file_path, ind):
    """Parse covariates and filter for desired sample IDs."""
    cov = pd.read_csv(cov_file_path, delimiter='\t', index_col="nwd_id").T
    cov = cov[ind]
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

def maf_filter(geno, ind, maf_thresh, ma_samp_thresh):
    """Exclude genotypes that do not meet minimum maf and minor allele sample
       thresholds."""
    # Create separate DataFrame to store allele counts 
    alleles = geno[ind] 
    
    # Record number of individuals with alternate allele for each SNP
    geno["alt_samp"] = alleles[alleles != 0].count(axis=1)
    
    # Record alternate allele frequency for each SNP
    n_ind = len(ind)
    geno["alt_af"] = alleles.sum(axis=1).divide(2 * n_ind)

    # Record minor allele frequency and count
    geno["ma_samp"] = pd.DataFrame([geno["alt_samp"], 
                                   n_ind - geno["alt_samp"]]).T.min(axis=1)
    geno["maf"] = pd.DataFrame([geno["alt_af"], 1 - geno["alt_af"]]).T.min(axis=1)

    # Remove rows that do not meet both minor allele thresholds
    geno = geno[(geno.ma_samp >= ma_samp_thresh) & 
                (geno.maf >= maf_thresh)][list(ind) + geno_info_cols]
    return(geno)

def perform_regression(geno, exp, cov):
    """Regress genotypes onto expression data."""
    def marginal_test(row, exp, cov):
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
    res[regression_stat_cols] = geno.apply(lambda x: 
        marginal_test(x, exp, cov), axis=1, result_type='expand')
    return(res)

def perm_adjust_pval(n, min_pval, geno, exp, ind):
    """Adjust p-values using permutations of genotype/expression data. Not
       currently used in pipeline."""
    def permutation(geno, exp, ind):
        # Randomize expression data
        randomized_ind = np.random.permutation(ind)
        exp = exp[randomized_ind]

        # Create dummy column names so that X and Y variables match for sm.OLS
        n_ind = len(ind)
        dummy_ind = list(range(n_ind))
        exp.columns, geno.columns = dummy_ind, dummy_ind 

        # Regress randomized expression on genotypes
        perm_min_pval = perform_regression(geno, exp, None)["pval"].min()
        return(perm_min_pval)

    null_pval_dist = pd.Series([permutation(geno, exp, ind) for _ in range(n)])
    adj_pval = (null_pval_dist[null_pval_dist < min_pval].count().sum() + 1) / (n + 1)
    return(adj_pval)

def remove_nan(geno, exp):
    """Remove individuals from regression if we do not have expression
       measurements for them."""
    all_ind = exp.columns
    exp_dropna = exp.dropna(axis=1)
    ind_dropna = exp_dropna.columns
    dropped_ind = all_ind.difference(ind_dropna)
    geno_dropna = geno.drop(dropped_ind, axis=1)
    return(geno_dropna, exp_dropna, ind_dropna)

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

    # If regression fails, set command to create blank output file
    fail_command = "touch {0}".format(args.out)

    # Record number of permutations to run
    n_perm = int(args.n_perm)
    
    # Load expression and genotype data for ALL individuals, then select the
    # individuals for the current regression (ind_IDs)
    ind_IDs = get_ind(args.samples)
    expression = parse_expression(args.phenotypes, ind_IDs)
    genotypes = parse_genotypes(args.genotypes, ind_IDs)
    if args.covariates is not None:
        covariates = parse_covariates(args.covariates, ind_IDs)
    else:
        covariates = None

    # Regress expression data on genotype data 
    if genotypes is None: 
        # Exit if there are no SNPs within 100 kb of the gene's TSS,
        # creating blank output file so that Snakemake doesn't throw a fit
        os.system(fail_command)
    else: 
        genotypes, expression, ind_IDs = remove_nan(genotypes, expression)
        genotypes = maf_filter(genotypes, ind_IDs, args.maf, args.ma_samples)
        if genotypes.empty:
            # Exit if there are no SNPs in the desired MAF window,
            # creating blank output file so that Snakemake doesn't throw a fit
            os.system(fail_command)
        else:
            results = perform_regression(genotypes, expression, covariates)
            results[geno_info_cols] = genotypes[geno_info_cols]

    # Adjust p-values with permutations if desired and write results to output
     if n_perm == 0:
         results.to_csv(args.out, sep='\t')
     else:
         # Select SNP with min pval (randomly choose in case of tie)
         candidate = results[results.pval == results.pval.min()].sample(1) 
         cand_pval = candidate["pval"].iloc[0]

         # If candidate SNP is significant, perform permutations
         if cand_pval < .10:
             perm_pval = perm_adjust_pval(n_perm, cand_pval, genotypes, expression, ind_IDs)
         else:
             candidate[["perm_pval", "n_perm"]] = [None, 0]
             candidate.to_csv(args.out, sep='\t')
         candidate["perm_pval"] = perm_pval
         candidate["n_perm"] = n_perm
         candidate.to_csv(args.out, sep='\t')

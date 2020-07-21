import pandas as pd
import argparse
import numpy as np
from call_eqtls import parse_genotypes
from map_ancestry_genotypes import find_match, anc_geno_string

def combine_pheno_files(file_list, hits):
    """Merges phenotypes for all significant eQTLs with the gene name and variant ID."""
    exp = pd.DataFrame()
    for f in file_list:
        tmp = pd.read_csv(f, sep='\t')
        exp = pd.concat([exp, tmp])
    exp = exp.set_index("gene_id")
    exp.index = exp.index.rename("Gene")
    exp["ID"] = exp.apply(lambda x: hits[hits.Gene == x.name]["ID"].iloc[0], axis=1)
    exp = exp.set_index("ID", append=True)
    return(exp)

def add_intercept(df, name, exp, geno):
    """Calculates mean expression among individuals homozygous for reference allele."""
    hom_ref_mask = geno.where(geno == '0-0').replace(to_replace='0-0', value=1)
    masked_exp = exp * hom_ref_mask
    hom_ref_mean = masked_exp.mean(axis=1).reset_index().rename(columns={0: "intercept_" + name})
    df = pd.merge(df, hom_ref_mean)
    return(df)

def combine_geno_files(file_list, hits):
    """Merges genotypes for all significant eQTLs with the gene name and variant ID."""
    hits = hits[["Gene", "ID"]]
    geno = pd.DataFrame()
    for f in file_list:
        curr_gene = f.split('/')[-1][:-7]
        curr_hits = hits[hits.Gene == curr_gene]
        tmp = parse_genotypes(f, concat=False)
        tmp = pd.merge(tmp, curr_hits, left_index=True, right_on="ID", how='inner')
        geno = pd.concat([geno, tmp])
    geno = geno.set_index(["Gene", "ID"])
    return(geno)

def ancestry_phase_genotypes(geno_df, swap, tracts=None):
    def european_phaser(geno):
        if swap:
            int_geno = [2 - (int(g[0]) + int(g[2])) for g in geno]
            anc_geno = ['{0}-0'.format(str(g)) for g in int_geno]
        else:
            int_geno = [int(g[0]) + int(g[2]) for g in geno]
            anc_geno = ['{0}-0'.format(str(g)) for g in int_geno]
        return(pd.Series(anc_geno))
    def african_phaser(col, col_tracts):
        variants = col.keys()
        ancestries = {'A': [], 'B': []}
        for hap in ancestries.keys():
            anc_list = []
            for var in variants:
                chrom = int(var[1].split('_')[0][3:])
                pos = int(var[1].split('_')[1])
                curr_tracts = col_tracts[(col_tracts.Haplotype == hap) & (col_tracts.Chrom == chrom)]
                mapped_tract_start = find_match(pos, curr_tracts)
                if type(mapped_tract_start) == int:
                    mapped_anc = curr_tracts.loc[mapped_tract_start].Ancestry
                    anc_list.append(mapped_anc)
                else: # Matching tract not found
                    anc_list.append(None)
            ancestries[hap] = anc_list
        zip_anc_geno = [i for i in zip(col.values, [anc for anc in zip(*list(ancestries.values()))])] # Zip together local ancestry information and genotypes
        anc_geno = [anc_geno_string(i, swap) for i in zip_anc_geno]
        return(pd.Series(anc_geno))
    if tracts is None: # European-American genotypes can be ancestry-phased w/o lanc information
        phased_df = geno_df.apply(european_phaser)
    else:
        phased_df = geno_df.apply(lambda x: african_phaser(x, tracts[tracts.NWDID == x.name]), axis=0, result_type='expand')
    phased_df.columns = geno_df.columns
    phased_df.index = geno_df.index
    return(phased_df)

def reshape_df(df, col_name):
    df = df.stack().reset_index().rename(columns={"level_2": "NWDID", 0: col_name})
    return(df)

def prepare_merged_df(afr_hit_path, eur_hit_path, tract_path, n_genes, swap_alleles, cov_paths):
    """Parses local ancestry tracts to determine ancestry-phased genotypes (i.e. # of alt. alleles
       in European ancestry region vs # of alt. alleles in African ancestry region). Resulting DataFrame
       has one row for each variant for each individual (i.e. total number of rows is nm, where 
       n is number of individuals and m is number of variants)."""
    # Merge significant eQTLs from regions of African ancestry in AA (also called Afr-hom) and regions of European ancestry in EA
    afr_hits = pd.read_csv(afr_hit_path, sep='\t')
    eur_hits = pd.read_csv(eur_hit_path, sep='\t')
    hits = pd.merge(afr_hits, eur_hits, on=["Gene", "ID"], suffixes = ("_Afr", "_Eur"))
    all_genes = list(hits["Gene"])

    # Read in local ancestry tracts (format: one row for each tract in each individual)
    tracts = pd.read_csv(tract_path, names=["Chrom", "Start", "Stop", "ID"], sep='\t')
    tracts[["NWDID", "Haplotype", "Ancestry"]] = tracts.apply(lambda x: pd.Series(x.ID.split('_')), axis=1) # Split tract ID into separate columns

    # Specify which genotype and phenotype files need to be read
    afr_geno_files = ["data/fastqtl_geno_input/Afr/" + g + ".vcf.gz" for g in all_genes]
    afr_pheno_files = ["data/fastqtl_pheno_input/Afr/" + g + ".txt" for g in all_genes]
    eur_geno_files = ["data/fastqtl_geno_input/Eur/" + g + ".vcf.gz" for g in all_genes]
    eur_pheno_files = ["data/fastqtl_pheno_input/Eur/" + g + ".txt" for g in all_genes]

    # Subsamples variants/genes/phenotypes based on user-specified argument, mostly for debugging purposes
    if n_genes is not None:
        n_genes = int(n_genes)
        afr_geno_files = afr_geno_files[:n_genes]
        afr_pheno_files = afr_pheno_files[:n_genes]
        eur_geno_files = eur_geno_files[:n_genes]
        eur_pheno_files = eur_pheno_files[:n_genes]

    # Reads in genotype files and creates m x n matrix of ancestry-phased genotypes (i.e. one row for each variant, 
    # one column for each individual
    afr_genotypes = combine_geno_files(afr_geno_files, hits)
    afr_genotypes = ancestry_phase_genotypes(afr_genotypes, swap_alleles, tracts)

    # Reads in phenotype files and creates m x n matrix of phenotypes (i.e. one row for each phenotype, 
    # one column for each individual
    afr_expression = combine_pheno_files(afr_pheno_files, hits)
    afr_expression = afr_expression[afr_genotypes.columns]

    # Reads in genotype files and creates m x n matrix of ancestry-phased genotypes (i.e. one row for each variant, 
    # one column for each individual
    eur_genotypes = combine_geno_files(eur_geno_files, hits)
    eur_genotypes = ancestry_phase_genotypes(eur_genotypes, swap_alleles)

    # Reads in phenotype files and creates m x n matrix of phenotypes (i.e. one row for each phenotype, 
    # one column for each individual
    eur_expression = combine_pheno_files(eur_pheno_files, hits)
    eur_expression = eur_expression[eur_genotypes.columns]

    # Adds population-specific intercept by calculating mean expression in reference homozygotes in each pop.
    # Now somewhat irrelevant/obsolete because we directly fit pop-specific intercepts (see below).
    hits = add_intercept(hits, "Afr", afr_expression, afr_genotypes)
    hits = add_intercept(hits, "Eur", eur_expression, eur_genotypes)

    # Reshapes genotype and phenotype dataframes such that there is one row for each variant/phenotype for each individual
    afr_genotypes = reshape_df(afr_genotypes, "Genotype")
    afr_genotypes["Race_AA"] = 1
    eur_genotypes = reshape_df(eur_genotypes, "Genotype")
    eur_genotypes["Race_AA"] = 0
    afr_expression = reshape_df(afr_expression, "Expression")
    eur_expression = reshape_df(eur_expression, "Expression")

    # Merge data and drop NaN entries
    expression = pd.concat([afr_expression, eur_expression])
    genotypes = pd.concat([afr_genotypes, eur_genotypes])
    merged_data = pd.merge(expression, genotypes)
    merged_data = pd.merge(merged_data, hits)

    # Merge covariate information if provided
    if cov_paths is not None:
        cov = pd.DataFrame()
        for path in cov_paths:
            tmp = pd.read_csv(path, sep='\t')
            cov = pd.concat([cov, tmp])
        merged_data = pd.merge(merged_data, cov)

    merged_data = merged_data.dropna(subset=["Genotype", "Expression", "effect_Afr", "effect_Eur", "Race_AA"])

    # Split ancestry-phased genotype into its constituent columns
    merged_data[["Genotype_Eur", "Genotype_Afr"]] = merged_data.apply(lambda x: pd.Series([int(x.Genotype[0]), int(x.Genotype[2])]), axis=1)

    return(merged_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tracts')
    parser.add_argument('--afr_hits')
    parser.add_argument('--eur_hits')
    parser.add_argument('--swap_ref_alt', action='store_true')
    parser.add_argument('--n_genes', default=None)
    parser.add_argument('--covariates', nargs='+')
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = prepare_merged_df(args.afr_hits, args.eur_hits, args.tracts, args.n_genes, args.swap_ref_alt, args.covariates)
    merged_data.to_csv(args.out, sep='\t', index=False)

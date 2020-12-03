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
    exp.index = exp.index.rename("gene")
    exp["ID"] = exp.apply(lambda x: hits[hits.gene == x.name]["ID"].iloc[0], axis=1)
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
    hits = hits[["gene", "ID"]]
    geno = pd.DataFrame()
    for f in file_list:
        curr_gene = f.split('/')[-1][:-7]
        curr_hits = hits[hits.gene == curr_gene]
        tmp = parse_genotypes(f, concat=False)
        tmp = pd.merge(tmp, curr_hits, left_index=True, right_on="ID", how='inner')
        geno = pd.concat([geno, tmp])
    geno = geno.set_index(["gene", "ID"])
    return(geno)

def map_local_ancestry(geno_df, tracts=None):
    """Returns a DataFrame with one column for every individual and one row for every SNP.
       Each entry corresponds to the comma-delimited local ancestry haplotypes for the 
       individual at the given SNP."""
    def european_mapper(geno):
        return(pd.Series(['CEU,CEU' for _ in geno]))
    def african_mapper(col, col_tracts):
        variants = col.keys() # Grab index, which consists of both gene and SNP ID
        ancestries = []
        for var in variants:
            curr_anc = []
            chrom = int(var[1].split('_')[0][3:])
            pos = int(var[1].split('_')[1])
            for hap in ['A', 'B']:
                curr_tracts = col_tracts[(col_tracts.haplotype == hap) & 
                                         (col_tracts.chrom == chrom)]
                mapped_tract_start = find_match(pos, curr_tracts)
                if type(mapped_tract_start) == int: # Matching tract successfully found
                    mapped_anc = curr_tracts.loc[mapped_tract_start].ancestry
                    curr_anc.append(mapped_anc)
                else: # Matching tract not found
                    curr_anc.append(None)
            anc_string = ','.join([str(anc) for anc in curr_anc])
            ancestries.append(anc_string)
        return(pd.Series(ancestries))
    if tracts is None:
        lanc_df = geno_df.apply(european_mapper)
    else:
        lanc_df = geno_df.apply(lambda x: african_mapper(x, tracts[tracts.nwd_id == x.name]), 
                                  axis=0, result_type='expand')
    lanc_df.columns = geno_df.columns
    lanc_df.index = geno_df.index
    return(lanc_df)

def sum_local_ancestry(lanc_df):
    """Returns a DataFrame with one column for every individual and one row for every SNP.
       Each entry corresponds to the number of African ancestry tracts for the individual 
       at the given SNP, summing across both haplotypes."""
    def count(anc_string):
        num_afr = 0
        for anc in anc_string.split(','):
            if anc == "YRI":
                num_afr += 1
        return(num_afr)
    return(lanc_df.apply(lambda col: [count(x) for x in col]))

def ancestry_phase_genotypes(geno_df, swap, lanc_df):
    """Returns a DataFrame with one column for every individual and one row for every SNP.
       Each entry corresponds to the ancestry-phased genotype (i.e. the number of alternative
       alleles on a European background and the number of alternative alleles on an African
       background, separated by a hyphen) for the individual at the given SNP."""
    def phase_helper(geno_col, lanc_col):
        # Zip together local ancestry information and genotypes
        zip_anc_geno = [i for i in zip(geno_col.values, lanc_col.values)] 
        anc_geno = [anc_geno_string(i, swap) for i in zip_anc_geno]
        return(pd.Series(anc_geno))
    phased_df = geno_df.apply(lambda x: phase_helper(x, lanc_df.loc[:,x.name]))
    phased_df.columns = geno_df.columns
    phased_df.index = geno_df.index
    return(phased_df)

def reshape_df(df, col_name):
    """Reshapes DataFrame into long format."""
    df = df.stack().reset_index().rename(columns={"level_2": "nwd_id", 0: col_name})
    return(df)

def prepare_merged_df(afr_hit_path, eur_hit_path, tract_path, n_genes, swap_alleles, cov_path):
    """Parses local ancestry tracts to determine ancestry-phased genotypes (i.e. # of alt. 
       alleles in European ancestry region vs # of alt. alleles in African ancestry region). 
       Resulting DataFrame has one row for each variant for each individual (i.e. total 
       number of rows is nm, where n is number of individuals and m is number of variants)."""
    # Merge significant eQTLs from regions of African ancestry in AA (also called Afr-hom) and 
    # regions of European ancestry in EA
    afr_hits = pd.read_csv(afr_hit_path, sep='\t')
    eur_hits = pd.read_csv(eur_hit_path, sep='\t')
    hits = pd.merge(afr_hits, eur_hits, on=["gene", "ID"], suffixes = ("_Afr", "_Eur"))
    all_genes = list(hits["gene"])

    # Read in local ancestry tracts (format: one row for each tract in each individual) and 
    # split tract info ID into separate columns
    tracts = pd.read_csv(tract_path, names=["chrom", "start", "stop", "ID"], sep='\t')
    tracts[["nwd_id", "haplotype", "ancestry"]] = tracts.apply(lambda x: pd.Series(x.ID.split('_')), 
                                                               axis=1) 

    # Specify which genotype and phenotype files need to be read
    afr_geno_files = ["data/QTL_geno_input/Afr/" + g + ".vcf.gz" for g in all_genes]
    afr_pheno_files = ["data/QTL_pheno_input/Afr/" + g + ".txt" for g in all_genes]
    eur_geno_files = ["data/QTL_geno_input/Eur/" + g + ".vcf.gz" for g in all_genes]
    eur_pheno_files = ["data/QTL_pheno_input/Eur/" + g + ".txt" for g in all_genes]

    # Subsamples variants/genes/phenotypes based on user-specified argument. Used mostly for 
    # debugging purposes in conjunction with the swap_alleles parameter.
    if n_genes is not None:
        n_genes = int(n_genes)
        afr_geno_files = afr_geno_files[:n_genes]
        afr_pheno_files = afr_pheno_files[:n_genes]
        eur_geno_files = eur_geno_files[:n_genes]
        eur_pheno_files = eur_pheno_files[:n_genes]

    # Reads in genotype files and creates m x n matrix of ancestry-phased genotypes (i.e. one
    # row for each variant and one column for each individual)
    afr_genotypes = combine_geno_files(afr_geno_files, hits)
    afr_lanc_info = map_local_ancestry(afr_genotypes, tracts)
    afr_summed_lanc = sum_local_ancestry(afr_lanc_info)
    afr_genotypes = ancestry_phase_genotypes(afr_genotypes, swap_alleles, afr_lanc_info)

    # Reads in phenotype files and creates m x n matrix of phenotypes (i.e. one row for each 
    # phenotype and one column for each individual)
    afr_expression = combine_pheno_files(afr_pheno_files, hits)
    afr_expression = afr_expression[afr_genotypes.columns]

    # Reads in genotype files and creates m x n matrix of ancestry-phased genotypes (i.e. one 
    # row for each variant and one column for each individual)
    eur_genotypes = combine_geno_files(eur_geno_files, hits)
    eur_lanc_info = map_local_ancestry(eur_genotypes)
    eur_summed_lanc = sum_local_ancestry(eur_lanc_info)
    eur_genotypes = ancestry_phase_genotypes(eur_genotypes, swap_alleles, eur_lanc_info)

    # Reads in phenotype files and creates m x n matrix of phenotypes (i.e. one row for each 
    # phenotype and one column for each individual)
    eur_expression = combine_pheno_files(eur_pheno_files, hits)
    eur_expression = eur_expression[eur_genotypes.columns]

    # Adds population-specific intercept by calculating mean expression in reference 
    # homozygotes in each pop. In the likelihood model, this is used as an initialization for 
    # the population-specific intercepts in order to fit the first iteration of delta.
    hits = add_intercept(hits, "Afr", afr_expression, afr_genotypes)
    hits = add_intercept(hits, "Eur", eur_expression, eur_genotypes)

    # Reshapes genotype and phenotype dataframes such that there is one row for each variant/
    # phenotype for each individual
    afr_genotypes = reshape_df(afr_genotypes, "genotype")
    eur_genotypes = reshape_df(eur_genotypes, "genotype")
    afr_lanc = reshape_df(afr_summed_lanc, "local_ancestry")
    eur_lanc = reshape_df(eur_summed_lanc, "local_ancestry")
    afr_expression = reshape_df(afr_expression, "expression")
    eur_expression = reshape_df(eur_expression, "expression")

    # Merge data and drop NaN entries
    expression = pd.concat([afr_expression, eur_expression])
    genotypes = pd.concat([afr_genotypes, eur_genotypes])
    lanc = pd.concat([afr_lanc, eur_lanc])
    merged_data = pd.merge(expression, genotypes)
    merged_data = pd.merge(merged_data, lanc)
    merged_data = pd.merge(merged_data, hits)
    merged_data = merged_data.dropna(subset=["genotype", "expression", "effect_Afr", 
                                             "effect_Eur", "local_ancestry"])

    # Merge covariate information if provided
    if cov_path is not None:
        cov = pd.read_csv(cov_path, sep='\t')
        merged_data = pd.merge(merged_data, cov)

    # Split ancestry-phased genotype into its constituent columns
    merged_data[["genotype_Eur", "genotype_Afr"]] = merged_data.apply(lambda x: pd.Series([int(x.genotype[0]), 
                                                                                           int(x.genotype[2])]), axis=1)
    
    # Split race column
    merged_data["race_Afr"] = merged_data["race"]
    merged_data["race_Eur"] = merged_data.apply(lambda x: 1 - x.race_Afr, axis=1)

    return(merged_data)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tracts')
    parser.add_argument('--afr_hits')
    parser.add_argument('--eur_hits')
    parser.add_argument('--swap_ref_alt', action='store_true')
    parser.add_argument('--n_genes', default=None)
    parser.add_argument('--covariates')
    parser.add_argument('--out')
    args = parser.parse_args()

    merged_data = prepare_merged_df(args.afr_hits, args.eur_hits, args.tracts, 
                                    args.n_genes, args.swap_ref_alt, args.covariates)
    merged_data.to_csv(args.out, sep='\t', index=False)

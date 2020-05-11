import sys
import os
import pandas as pd
import argparse
import gzip
import numpy as np
import statsmodels.api as sm
from call_eqtls import parse_genotypes
from map_ancestry_genotypes import find_match, anc_geno_string

def combine_pheno_files(file_list, hits):
    exp = pd.DataFrame()
    for f in file_list:
        tmp = pd.read_csv(f, sep='\t')
        exp = pd.concat([exp, tmp])
    exp = exp.set_index("gene_id")
    exp.index = exp.index.rename("Gene")
    exp["ID"] = exp.apply(lambda x: hits[hits.Gene == x.name]["ID"].iloc[0], axis=1)
    exp = exp.set_index("ID", append=True)
    return(exp)

def mean_center_expression(exp, geno):
    hom_ref_mask = geno.where(geno == '0-0').replace(to_replace='0-0', value=1)
    masked_exp = exp * hom_ref_mask
    hom_ref_mean = masked_exp.mean(axis=1)
    centered_exp = exp.subtract(hom_ref_mean, axis=0)
    return(centered_exp)

def combine_geno_files(file_list, hits):
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

def ancestry_phase_genotypes(geno_df, tracts=None):
    def european_phaser(geno):
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
        anc_geno = [anc_geno_string(i) for i in zip_anc_geno]
        return(pd.Series(anc_geno))
    if tracts is None: # European-American genotypes can be ancestry-phased w/o lanc information
        phased_df = geno_df.apply(european_phaser)
    else:
        phased_df = geno_df.apply(lambda x: african_phaser(x, tracts[tracts.NWDID == x.name]), axis=0, result_type='expand')
    phased_df.columns = geno_df.columns
    phased_df.index = geno_df.index
    return(phased_df)        

def expand_genotypes(row):
    n_ind = length(row.values)
    expanded = np.zeros((n_ind, n_ind *2))
    i = 0
    for geno in row.values:
        expanded[i * 2] = int(geno[0])
        expanded[i * 2 + 1] = int(geno[2])
    expanded = pd.DataFrame(expanded)
    expanded.columns = [k + anc for anc in ['CEU', 'YRI'] for k in row.keys()]
    expanded.index = [row.name] * n_ind
    return(expanded)

def likelihood(df, delta):
    def rowwise_likelihood(row, delta):
        exp = row.Expression
        additive_effects = row.Genotype_Afr * row.effect_Afr + row.Genotype_Eur * row.effect_Eur
        interaction = delta * (row.effect_Afr - row.effect_Eur) * row.Genotype_Eur * row.Race_AA
        row_likelihood = (exp - additive_effects - interaction) ** 2
        return(row_likelihood)
    df_likelihood = df.apply(lambda x: rowwise_likelihood(x, delta), axis=1)
    return(df_likelihood.sum())
        
def reshape_df(df, col_name):
    df = df.stack().reset_index().rename(columns={"level_2": "NWDID", 0: col_name})
    return(df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--tracts')
    parser.add_argument('--afr_hits')
    parser.add_argument('--eur_hits')
    parser.add_argument('--out')
    args = parser.parse_args()
    
    afr_hits = pd.read_csv(args.afr_hits, sep='\t')
    eur_hits = pd.read_csv(args.eur_hits, sep='\t')
    hits = pd.merge(afr_hits, eur_hits, on=["Gene", "ID"], suffixes = ("_Afr", "_Eur"))
    all_genes = list(hits["Gene"])

    tracts = pd.read_csv(args.tracts, names=["Chrom", "Start", "Stop", "ID"], sep='\t')
    tracts[["NWDID", "Haplotype", "Ancestry"]] = tracts.apply(lambda x: pd.Series(x.ID.split('_')), axis=1) # Split tract ID into separate columns
    
    afr_geno_files = ["data/fastqtl_geno_input/Afr/" + g + ".vcf.gz" for g in all_genes]
    afr_pheno_files = ["data/fastqtl_pheno_input/Afr/" + g + ".txt" for g in all_genes]
    eur_geno_files = ["data/fastqtl_geno_input/Eur/" + g + ".vcf.gz" for g in all_genes]
    eur_pheno_files = ["data/fastqtl_pheno_input/Eur/" + g + ".txt" for g in all_genes]

    afr_genotypes = combine_geno_files(afr_geno_files, hits)
    afr_genotypes = ancestry_phase_genotypes(afr_genotypes, tracts)
    # expanded_afr_genotypes.index = [i for i in afr_genotypes.index for _ in afr_genotypes.shape[1]]

    afr_expression = combine_pheno_files(afr_pheno_files, hits)
    afr_expression = mean_center_expression(afr_expression, afr_genotypes) 
    afr_expression = afr_expression[afr_genotypes.columns]
    # expanded_afr_expression.index = [i for i in afr_expression.index for _ in afr_expression.shape[1]]

    eur_genotypes = combine_geno_files(eur_geno_files, hits)
    eur_genotypes = ancestry_phase_genotypes(eur_genotypes)
    # expanded_eur_genotypes.index = [i for i in eur_genotypes.index for _ in eur_genotypes.shape[1]]

    eur_expression = combine_pheno_files(eur_pheno_files, hits)
    eur_expression = mean_center_expression(eur_expression, eur_genotypes) 
    eur_expression = eur_expression[eur_genotypes.columns]
    # expanded_eur_expression.index = [i for i in eur_expression.index for _ in eur_expression.shape[1]]

    afr_genotypes = reshape_df(afr_genotypes, "Genotype")
    afr_genotypes["Race_AA"] = 1
    eur_genotypes = reshape_df(eur_genotypes, "Genotype")
    eur_genotypes["Race_AA"] = 0 

    afr_expression = reshape_df(afr_expression, "Expression")
    eur_expression = reshape_df(eur_expression, "Expression")

    expression = pd.concat([afr_expression, eur_expression])
    genotypes = pd.concat([afr_genotypes, eur_genotypes])
    merged_data = pd.merge(expression, genotypes)
    merged_data = pd.merge(merged_data, hits)
    merged_data = merged_data.dropna(subset=["Genotype", "Expression", "effect_Afr", "effect_Eur", "Race_AA"])
    merged_data[["Genotype_Eur", "Genotype_Afr"]] = merged_data.apply(lambda x: pd.Series([int(x.Genotype[0]), int(x.Genotype[2])]), axis=1)
    # genotypes = pd.concat([expand_genotypes(row) for i, row in genotypes.iterrows()])

    delta_list = [i / 100 for i in range(0, 101)]
    likelihood_list = []
    for d in delta_list:
        likelihood_list.append(likelihood(merged_data, d))
    res = pd.DataFrame({"Delta": delta_list, "Likelihood": likelihood_list})
    res.to_csv(args.out, sep='\t')

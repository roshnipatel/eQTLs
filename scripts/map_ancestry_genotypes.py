import pandas as pd
import numpy as np
import argparse
from call_eqtls import parse_genotypes

def find_match(x, df):
    for idx, row in df.iterrows():
        if row.Stop > x:
            if row.Start < x:
                return(idx)
            return('before')
    return('after')

def anc_geno_string(tup):
    geno = {'CEU': 0, 'YRI': 0}
    for i in range(2):
        anc = tup[1][i]
        if anc != 0:
            count = int(tup[0].split('|')[i])
            geno[anc] += count
        else:
            return(None)
    geno_str = '-'.join([str(i) for i in [geno['CEU'], geno['YRI']]])
    return(geno_str)

def find_SNP_ancestry(pos, tracts):
    n_snp = len(pos)
    anc = []
    first_idx = find_match(pos[0], tracts) # Find tract corresponding to first SNP position
    if first_idx == 'after': # All SNPs fall after the last known ancestry tract
        anc = [0] * n_snp
    else:
        last_idx = find_match(pos[n_snp - 1], tracts)
        if (first_idx == 'before') & (last_idx == 'before'): # All SNPs fall in a gap between two ancestry tracts
            anc = [0] * n_snp
        elif (first_idx == 'before') & (last_idx == 'after'): # Something is wrong
            print("Houston, we have a problem.")
        elif first_idx == last_idx: # All SNPs fall in the same ancestry tract
            anc = [tracts.loc[first_idx].Ancestry] * n_snp
        else: # Not all SNPs fall in the same ancestry tract; iterate 
            for p in pos:
                curr_idx = find_match(p, tracts)
                anc.append(curr_idx)
    return(anc)

def map_ancestry(row, tracts):
    tracts = tracts[tracts.NWDID == row.name]
    positions = [int(k.split('_')[1]) for k in row.keys()]
    ancestries = {'A': [], 'B': []}
    for hap in ancestries.keys():
        curr_tracts = tracts[tracts.Haplotype == hap]
        ancestries[hap] = find_SNP_ancestry(positions, curr_tracts)
    zip_lanc_geno = [i for i in zip(row.values, [anc for anc in zip(*list(ancestries.values()))])] # Zip together local ancestry information and genotypes
    lanc_geno = [anc_geno_string(i) for i in zip_lanc_geno]
    return(pd.Series(lanc_geno))

def split_tract_ID(row):
    return(pd.Series(row.ID.split('_')))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--lanc', help='bed file of local ancestry tracts')
    parser.add_argument('--genotypes', help='VCF of genotypes to map onto tracts. All SNPs should be on same chromosome.')
    parser.add_argument('--out', help='output file')
    args = parser.parse_args()
    
    genotypes = parse_genotypes(args.genotypes, concat=False).T

    tracts = pd.read_csv(args.lanc, names=["Chrom", "Start", "Stop", "ID"], sep='\t')
    tracts = tracts[tracts.Chrom == int(genotypes.columns[0].split('_')[0][3:])] # Subset for the current chromosome
    tracts[["NWDID", "Haplotype", "Ancestry"]] = tracts.apply(split_tract_ID, axis=1) # Split tract ID into separate columns
    
    lanc_cols = [col + "_CEU-YRI" for col in genotypes.columns]
    genotypes[lanc_cols] = genotypes.apply(lambda x: map_ancestry(x, tracts), axis=1)
    genotypes.to_csv(args.out, sep='\t')

import pandas as pd
import gzip
import argparse

def parse_genotypes(vcf_path):
    """Lightly adapted from function in estimate_effect_sizes.py to store
       the information we usually don't care about."""
    header = []
    geno_list = []
    with gzip.open(vcf_path, 'rt') as geno_file:
        for line in geno_file:
            if line[:2] == "##":
                header.append(line.strip())
                continue
            if line[0] == "#":
                line = line.split()
                all_cols = line
                ind_cols = line[9:]
            else:
                line = line.split()
                geno_list.append(line)
    geno_df = pd.DataFrame(geno_list, columns=all_cols)
    return(header, geno_df, ind_cols)

def drop_genotypes(geno, keep=None):
    if not keep:
        return(".|.")
    elif keep == "A":
        return("{0}|.".format(geno[0]))
    elif keep == "B":
        return(".|{0}".format(geno[2]))

def write_vcf(header, geno_df, filepath): 
    header_string = '\n'.join(header)
    df_string = geno_df.to_csv(sep='\t', index=False)
    with gzip.open(filepath, 'wt') as f:
        f.write(header_string)
        f.write('\n')
        f.write(df_string)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--tracts')
    parser.add_argument('--anc')
    parser.add_argument('--gene')
    parser.add_argument('--out')
    args = parser.parse_args()

    tracts = pd.read_csv(args.tracts, delimiter='\t',
                        names=["chrom", "anc_start", "anc_stop", "info", "gene_chrom",
                               "gene_start", "gene_stop", "gene_id", "overlap"],
                        dtype={"chrom": int, "anc_start": int, "anc_stop": int,
                               "info": str, "gene_chrom": str, "gene_start": int,
                               "gene_stop": int, "gene_id": str, "overlap": int})
    tracts = tracts.loc[tracts.gene_id == args.gene]
    tracts["nwd_id"] = tracts.apply(lambda row: row.info.split('_')[0], axis=1)
    tracts["hap"] = tracts.apply(lambda row: row.info.split('_')[1], axis=1)
    tracts["anc"] = tracts.apply(lambda row: row.info.split('_')[2], axis=1)
    
    # Filter for individuals with 200 kb of desired ancestry
    if args.anc == "Afr":
        curr_anc = "YRI"
    elif args.anc == "Eur":
        curr_anc = "CEU"
    tracts_to_keep = tracts.loc[(tracts.anc == curr_anc) & (tracts.overlap == 200000)]

    vcf_header, genotypes, geno_ind = parse_genotypes(args.vcf)
    for ind in geno_ind:
        curr_tracts = tracts_to_keep.loc[tracts_to_keep.nwd_id == ind]
        if curr_tracts.empty:
            genotypes.loc[:,ind] = genotypes.loc[:,ind].apply(drop_genotypes)
        elif curr_tracts.shape[0] == 1:
            hap_to_keep = curr_tracts.hap.iloc[0]
            genotypes.loc[:,ind] = genotypes.loc[:,ind].apply(lambda x: drop_genotypes(x, keep=hap_to_keep))

    write_vcf(vcf_header, genotypes, args.out)

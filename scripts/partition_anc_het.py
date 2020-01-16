import argparse
import sys
import pandas as pd
from combine_tracts import parse_tract_files 

parser = argparse.ArgumentParser()
parser.add_argument('--hits')
parser.add_argument('--samples')
parser.add_argument('--tract_dir')
parser.add_argument('--chrom_lengths')
parser.add_argument('--sample_out')
parser.add_argument('--region_out')
args = parser.parse_args()

# Load dataframe of significant hits
hits = pd.read_csv(args.hits, delimiter='\t')
SNP_list = hits["variant"]

# Load dataframe of samples
samples = pd.read_csv(args.samples, delimiter='\t')

# Load dataframe of ancestry tracts
tracts_df = parse_tract_files(args.tract_dir)
tracts_df["NWDID"] = tracts_df["Info"].apply(lambda x: x.split('_')[0])
tracts_df = pd.merge(tracts_df, samples, how='right')[["NWDID", "Chrom", "Start", "Stop", "Anc", "Info"]]

# Load dataframe of chromosome lengths
chrom_lengths = pd.read_csv(args.chrom_lengths, delimiter='\t', names=["Chrom", "Length"])

# Iterate through VCF.
for _, SNP in SNP_list.iteritems(): 
    SNP_split = SNP.split('_')
    chrom = int(SNP_split[0][3:])
    pos = int(SNP_split[1])
    
    # Identify individuals where SNP falls in region of European ancestry on one 
    # haplotype, and a region of African ancestry on the other haplotype. Exclude 
    # SNPs that are within 200kb of ancestry breakpoint (required due to the 
    # specific parameters of partition_samples.py).
    start_SNP_window = max(0, pos - 200000)
    end_SNP_window = min(pos + 200000, chrom_lengths[chrom_lengths["Chrom"] == str(chrom)]["Length"].iloc[0])
    eur_tracts = tracts_df[(tracts_df.Chrom == chrom) & (tracts_df.Anc == 'CEU') & (tracts_df.Start <= start_SNP_window) & (tracts_df.Stop >= end_SNP_window)]
    afr_tracts = tracts_df[(tracts_df.Chrom == chrom) & (tracts_df.Anc == 'YRI') & (tracts_df.Start <= start_SNP_window) & (tracts_df.Stop >= end_SNP_window)]
    anc_het_ind = pd.merge(afr_tracts, eur_tracts, how='inner', on="NWDID")
    
    # Write individuals to file.
    anc_het_ind.to_csv(args.sample_out + SNP + ".txt", columns=["NWDID"], header=False, index=False) 
    
    # Write the 2 Mb region around SNP to BED-formatted file.
    region_size = 1000000
    start_region = max(0, pos - region_size)
    end_region = min(pos + region_size, chrom_lengths[chrom_lengths["Chrom"] == str(chrom)]["Length"].iloc[0])
    f = open(args.region_out + SNP + ".txt", 'w')
    f.write("chr{0}\t{1}\t{2}\t{3}".format(chrom, start_region, end_region, SNP))
    f.close()

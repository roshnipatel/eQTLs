import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--intersect')
parser.add_argument('--samples')
parser.add_argument('--metadata')
parser.add_argument('--genes')
parser.add_argument('--out')
args = parser.parse_args()

cis_window = 2000000
overlap = 1800000

samples = pd.read_csv(args.samples, delimiter='\t')
genes = pd.read_csv(args.genes, delimiter='\t', names=["Chrom", "Start", "Stop", "GeneID"])["GeneID"]

# Filter tracts data for individuals and genes in our sample dataset.
tracts = pd.read_csv(args.intersect, delimiter='\t',
                    names=["Chrom", "AncStart", "AncStop", "Ind_Hapl", "GeneChrom", "GeneStart", "GeneStop", "GeneID", "Overlap"],
                    dtype={"Chrom": int, "AncStart": int, "AncStop": int, "Ind_Hapl": str, "GeneChrom": str, "GeneStart": int, "GeneStop": int, "GeneID": str, "Overlap": int})
tracts = pd.merge(genes, tracts, how='inner')
tracts["NWDID"] = tracts.apply(lambda row: row.Ind_Hapl[:-2], axis=1)
tracts = pd.merge(samples, tracts, how='left')[["GeneID", "NWDID", "Overlap"]]

# Identify Afr-Am individuals that have African ancestry overlapping gene's cis
# window for both chromosomes. Filter for individuals in our sample dataset.
Afr_tracts = tracts[tracts["Overlap"] > overlap]
Afr_tracts = Afr_tracts.groupby(["GeneID", "NWDID"]).size()
Afr_tracts = Afr_tracts[Afr_tracts > 1]
Afr_tracts = Afr_tracts.reset_index(name='counts')
Afr_tracts = Afr_tracts.drop(['counts'], axis=1)

# Write Afr-Am sample IDs to file and store sample sizes for each gene
gene_counts = {}
for gene, df in Afr_tracts.groupby("GeneID"):
    gene_counts[gene] = df.shape[0]
    df.to_csv(args.out + "/estimation/Afr/" + gene + ".txt", header=False, index=False, columns=["NWDID"])  

# Partition Eur individuals into ascertainment and estimation set, such that size of 
# Eur estimation set matches size of Afr-Am estimation set.
metadata = pd.read_csv(args.metadata)[["NWDID", "race1c"]]
metadata = metadata.drop_duplicates()
merged = pd.merge(samples, metadata, how='left')
Eur_IDs = merged[merged.race1c == 1]["NWDID"]

for gene, count in gene_counts.items():
    est = Eur_IDs.sample(n = count)
    asc = Eur_IDs.drop(est.index)
    est.to_csv(args.out + "/estimation/Eur/" + gene + ".txt", header=False, index=False)
    asc.to_csv(args.out + "/ascertainment/Eur/" + gene + ".txt", header=False, index=False)

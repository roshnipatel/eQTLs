import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--intersect')
parser.add_argument('--samples')
parser.add_argument('--metadata')
parser.add_argument('--genes')
parser.add_argument('--out')
args = parser.parse_args()

cis_window = 200000

samples = pd.read_csv(args.samples, delimiter='\t')
genes = pd.read_csv(args.genes, delimiter='\t', names=["Chrom", "Start", "Stop", "GeneID"])["GeneID"]

# Filter tracts data for individuals and genes in our sample dataset.
tracts = pd.read_csv(args.intersect, delimiter='\t',
                    names=["Chrom", "AncStart", "AncStop", "Info", "GeneChrom", "GeneStart", "GeneStop", "GeneID", "Overlap"],
                    dtype={"Chrom": int, "AncStart": int, "AncStop": int, "Info": str, "GeneChrom": str, "GeneStart": int, "GeneStop": int, "GeneID": str, "Overlap": int})
tracts = pd.merge(genes, tracts, how='inner')
tracts["NWDID"] = tracts.apply(lambda row: row.Info.split('_')[0], axis=1)
tracts["Anc"] = tracts.apply(lambda row: row.Info.split('_')[2], axis=1)
tracts = pd.merge(samples, tracts, how='left')[["GeneID", "NWDID", "Anc", "Overlap"]]

# Identify tracts of African ancestry
Afr_tracts = tracts[(tracts["Overlap"] == cis_window) & (tracts["Anc"] == 'YRI')]
Afr_tracts = Afr_tracts.drop(['Anc'], axis=1)
Afr_tracts = Afr_tracts.groupby(["GeneID", "NWDID"]).size()

# Identify tracts of European ancestry
Eur_tracts = tracts[(tracts["Overlap"] == cis_window) & (tracts["Anc"] == 'CEU')]
Eur_tracts = Eur_tracts.drop(["Anc"], axis=1)
Eur_tracts = Eur_tracts.groupby(["GeneID", "NWDID"]).size()

# Identify Afr-Am individuals that have African ancestry overlapping gene's cis
# window for both chromosomes.
hom_Afr_tracts = Afr_tracts[Afr_tracts > 1]
hom_Afr_tracts = hom_Afr_tracts.reset_index(name='counts')
hom_Afr_tracts = hom_Afr_tracts.drop(['counts'], axis=1)

# Write Afr-Am sample IDs to file and store sample sizes for each gene
gene_counts = {}
for gene, df in hom_Afr_tracts.groupby("GeneID"):
    gene_counts[gene] = df.shape[0]
    df.to_csv(args.out + "/estimation/Afr/" + gene + ".txt", header=False, index=False, columns=["NWDID"])  

# Identify Afr-Am individuals that have one chromosome with African ancestry
# overlapping gene's cis window and one chromosome with European ancestry 
# overlapping gene's cis window.
het_Afr_tracts = Afr_tracts[Afr_tracts == 1]
het_Afr_tracts = het_Afr_tracts.reset_index(name='counts')
het_Eur_tracts = Eur_tracts[Eur_tracts == 1]
het_Eur_tracts = het_Eur_tracts.reset_index(name='counts')
anc_het_tracts = pd.merge(het_Afr_tracts, het_Eur_tracts, how='inner')

# Write ancestry-heterozygous sample IDs to file
for gene, df in anc_het_tracts.groupby("GeneID"):
    df.to_csv(args.out + "/estimation/het/" + gene + ".txt", header=False, index=False, columns=["NWDID"])

tracts, het_Afr_tracts, het_Eur_tracts, hom_Afr_tracts, anc_het_tracts, Afr_tracts, Eur_tracts = None, None, None, None, None, None, None

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

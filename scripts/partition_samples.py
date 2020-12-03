import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--intersect')
parser.add_argument('--afr_samples')
parser.add_argument('--eur_samples')
parser.add_argument('--afr_validation_run', action='store_true')
parser.add_argument('--genes')
parser.add_argument('--out_dir')
parser.add_argument('--partition')
args = parser.parse_args()

cis_window = 200000

afr_samples = pd.read_csv(args.afr_samples, delimiter='\t')
eur_samples = pd.read_csv(args.eur_samples, delimiter='\t')
genes = pd.read_csv(args.genes, delimiter='\t', names=["chrom", "start", "stop", "gene_id"])["gene_id"]

# Filter tracts data for individuals and genes in our sample dataset.
tracts = pd.read_csv(args.intersect, delimiter='\t',
                    names=["chrom", "anc_start", "anc_stop", "info", "gene_chrom", "gene_start", "gene_stop", "gene_id", "overlap"],
                    dtype={"chrom": int, "anc_start": int, "anc_stop": int, "info": str, "gene_chrom": str, "gene_start": int, "gene_stop": int, "gene_id": str, "overlap": int})
tracts = pd.merge(genes, tracts, how='inner')
tracts["nwd_id"] = tracts.apply(lambda row: row.info.split('_')[0], axis=1)
tracts["anc"] = tracts.apply(lambda row: row.info.split('_')[2], axis=1)
tracts = pd.merge(afr_samples, tracts, how='left')[["gene_id", "nwd_id", "anc", "overlap"]]

# Identify tracts of African ancestry
Afr_tracts = tracts[(tracts["overlap"] == cis_window) & (tracts["anc"] == 'YRI')]
Afr_tracts = Afr_tracts.drop(['anc'], axis=1)
Afr_tracts = Afr_tracts.groupby(["gene_id", "nwd_id"]).size()

# Identify tracts of European ancestry
Eur_tracts = tracts[(tracts["overlap"] == cis_window) & (tracts["anc"] == 'CEU')]
Eur_tracts = Eur_tracts.drop(["anc"], axis=1)
Eur_tracts = Eur_tracts.groupby(["gene_id", "nwd_id"]).size()

# Identify Afr-Am individuals that have African ancestry overlapping gene's cis
# window for both chromosomes.
hom_Afr_tracts = Afr_tracts[Afr_tracts > 1]
hom_Afr_tracts = hom_Afr_tracts.reset_index(name='counts')
hom_Afr_tracts = hom_Afr_tracts.drop(['counts'], axis=1)

if args.afr_validation_run:
    for gene, df in hom_Afr_tracts.groupby("gene_id"):
        df.sample(n = 50).to_csv(args.out_dir + "/reestimation_validation/Afr/" + gene + ".txt", header=False, index=False, columns=["nwd_id"])

# Write Afr-Am sample IDs to file and store sample sizes for each gene
gene_counts = {}
for gene, df in hom_Afr_tracts.groupby("gene_id"):
    gene_counts[gene] = df.shape[0]
    df.to_csv(args.out_dir + "/reestimation_primary/Afr/" + gene + ".txt", header=False, index=False, columns=["nwd_id"])  

# Identify Afr-Am individuals that have one chromosome with African ancestry
# overlapping gene's cis window and one chromosome with European ancestry 
# overlapping gene's cis window.
het_Afr_tracts = Afr_tracts[Afr_tracts == 1]
het_Afr_tracts = het_Afr_tracts.reset_index(name='counts')
het_Eur_tracts = Eur_tracts[Eur_tracts == 1]
het_Eur_tracts = het_Eur_tracts.reset_index(name='counts')
anc_het_tracts = pd.merge(het_Afr_tracts, het_Eur_tracts, how='inner')

# Write ancestry-heterozygous sample IDs to file
for gene, df in anc_het_tracts.groupby("gene_id"):
    df.to_csv(args.out_dir + "/reestimation_primary/het/" + gene + ".txt", header=False, index=False, columns=["nwd_id"])

# Initialize matrix to store ascertainment information
all_ind = pd.concat([eur_samples.nwd_id, afr_samples.nwd_id])
all_genes = gene_counts.keys()
partition_matrix = pd.DataFrame(np.zeros((len(all_genes), len(all_ind))))
partition_matrix.columns = all_ind
partition_matrix.index = all_genes

# Partition Eur individuals into ascertainment and estimation set, such that size of 
# Eur estimation set matches size of Afr-Am estimation set.
Eur_IDs = eur_samples["nwd_id"]
for gene, count in gene_counts.items():
    est = Eur_IDs.sample(n = count)
    rest = Eur_IDs.drop(est.index)
    val = rest.sample(n = 50)
    asc = rest.drop(val.index)
    partition_matrix.loc[gene, asc] = 1
    val.to_csv(args.out_dir + "/reestimation_validation/Eur/" + gene + ".txt", header=False, index=False)
    est.to_csv(args.out_dir + "/reestimation_primary/Eur/" + gene + ".txt", header=False, index=False)
    asc.to_csv(args.out_dir + "/ascertainment/Eur/" + gene + ".txt", header=False, index=False)
partition_matrix.to_csv(args.partition, sep='\t', header=True, index=True)

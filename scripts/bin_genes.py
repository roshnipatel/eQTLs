# Run in py36 conda environment
import pandas as pd

gene_info = pd.read_csv("data/hakha_annotations.txt", sep='\t')
gene_ensg = pd.read_csv("data/gene_ENSG_table.txt",sep='\t').loc[:,["GeneSymbol", "gene"]]
gene_info = pd.merge(gene_info, gene_ensg, on="gene")

# TF bins
gene_info.loc[gene_info.TF == 1].to_csv("data/bins/gene_bin_TF_1.txt", columns=["GeneSymbol"])
gene_info.loc[gene_info.TF == 0].to_csv("data/bins/gene_bin_TF_0.txt", columns=["GeneSymbol"])

# Connectedness bins
gene_info.loc[gene_info.connect_rank == 1].to_csv("data/bins/gene_bin_connectedness_high.txt", columns=["GeneSymbol"])
gene_info.loc[gene_info.connect_rank == 11].to_csv("data/bins/gene_bin_connectedness_low.txt", columns=["GeneSymbol"])
gene_info.loc[(gene_info.connect_rank > 1) & (gene_info.connect_rank < 11)].to_csv("data/bins/gene_bin_connectedness_med.txt", columns=["GeneSymbol"])

# pLI bins
pLI_quantiles = gene_info.pLI.describe()
gene_info.loc[gene_info.pLI < pLI_quantiles["25%"]].to_csv("data/bins/gene_bin_pLI_q1.txt", columns=["GeneSymbol"])
gene_info.loc[(gene_info.pLI >= pLI_quantiles["25%"]) & (gene_info.pLI < pLI_quantiles["50%"])].to_csv("data/bins/gene_bin_pLI_q2.txt", columns=["GeneSymbol"])
gene_info.loc[(gene_info.pLI >= pLI_quantiles["50%"]) & (gene_info.pLI < pLI_quantiles["75%"])].to_csv("data/bins/gene_bin_pLI_q3.txt", columns=["GeneSymbol"])
gene_info.loc[gene_info.pLI >= pLI_quantiles["75%"]].to_csv("data/bins/gene_bin_pLI_q4.txt", columns=["GeneSymbol"])

# EDS bins
EDS_quantiles = gene_info.EDS.describe()
gene_info.loc[gene_info.EDS < EDS_quantiles["25%"]].to_csv("data/bins/gene_bin_EDS_q1.txt", columns=["GeneSymbol"])
gene_info.loc[(gene_info.EDS >= EDS_quantiles["25%"]) & (gene_info.EDS < EDS_quantiles["50%"])].to_csv("data/bins/gene_bin_EDS_q2.txt", columns=["GeneSymbol"])
gene_info.loc[(gene_info.EDS >= EDS_quantiles["50%"]) & (gene_info.EDS < EDS_quantiles["75%"])].to_csv("data/bins/gene_bin_EDS_q3.txt", columns=["GeneSymbol"])
gene_info.loc[gene_info.EDS >= EDS_quantiles["75%"]].to_csv("data/bins/gene_bin_EDS_q4.txt", columns=["GeneSymbol"])

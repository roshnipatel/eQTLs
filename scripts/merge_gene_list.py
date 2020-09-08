import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bed', nargs='+', help="bed files for which gene lists need to be merged")
parser.add_argument('--out', help="output file name")
args = parser.parse_args()

cols_to_keep = ["#chr", "start", "end", "gene_id"]

bed_df = pd.DataFrame()
for filepath in args.bed:
    tmp = pd.read_csv(filepath, sep='\t')[cols_to_keep]
    bed_df = pd.concat([bed_df, tmp])

shared_genes = bed_df.groupby("gene_id").size()[bed_df.groupby("gene_id").size() == 2].index
bed_df = bed_df[bed_df.gene_id.isin(shared_genes)].drop_duplicates()
bed_df = bed_df[bed_df["#chr"] != "chrX"]
bed_df.to_csv(args.out, index=False, sep='\t')

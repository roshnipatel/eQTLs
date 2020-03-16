import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--bed', nargs='+', help="bed files for which gene lists need to be merged")
parser.add_argument('--out', help="output file name")
args = parser.parse_args()

cols_to_keep = ["#chr", "start", "end", "gene_id"]

first = True
for filepath in args.bed:
    if first:
        first = False
        bed_df = pd.read_csv(filepath, sep='\t')[cols_to_keep]
    else:
        tmp = pd.read_csv(filepath, sep='\t')[cols_to_keep]
        bed_df = pd.merge(bed_df, tmp)

bed_df.to_csv(args.out, index=False, sep='\t')

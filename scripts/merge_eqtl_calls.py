import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--dir', help='Directory with individual eQTL files')
parser.add_argument('--out', help='Output file name')
args = parser.parse_args()

eqtl_files = os.listdir(args.dir)

merged = pd.DataFrame()
i = 0
for f in eqtl_files:
    try:
        tmp = pd.read_csv(args.dir + f, sep='\t')
    except pd.errors.EmptyDataError: # Catch error that results from an empty results file (e.g. if there were no SNPs to test)
        continue
    tmp["gene"] = f[:-4]
    merged = pd.concat([merged, tmp])
    print("Merged file number {0}".format(i))
    i += 1

merged.to_csv(args.out, sep='\t', index=False)

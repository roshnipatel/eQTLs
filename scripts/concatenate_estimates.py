import pandas as pd
import argparse
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--dir', help='Directory with individual gene files')
    parser.add_argument('--out', help='Output file name')
    args = parser.parse_args()
    
    eqtl_files = os.listdir(args.dir)
    
    # Concatenate estimates of eQTL/SNP effect sizes for all genes
    concatenated = pd.DataFrame()
    i = 0
    for f in eqtl_files:
        try:
            tmp = pd.read_csv(args.dir + f, sep='\t')
        except pd.errors.EmptyDataError: 
            # Catch error that results from an empty results file (e.g. no SNPs to test)
            continue
        tmp["gene"] = f[:-4]
        concatenated = pd.concat([merged, tmp])
        i += 1
    
    merged.to_csv(args.out, sep='\t', index=False)

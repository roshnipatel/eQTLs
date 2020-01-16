# Combines multiple local ancestry files (based on output from LocalAncestry
# pipeline) into a single BED-formatted file. Important: output is NOT sorted.

import argparse
import pandas as pd
import os
import numpy as np

def parse_tract_files(tract_dir):
    files = os.listdir(tract_dir)
    df = pd.DataFrame(columns=["Chrom", "Start", "Stop", "Anc", "Info"])
    for f in files:
        ind_ID = f.split('.')[2] # WILL NEED TO CHANGE THIS IF YOU UPDATE THE OUTPUT OF LOCAL ANCESTRY PIPELINE
        hapl = f.split('.')[3] # WILL NEED TO CHANGE THIS IF YOU UPDATE THE OUTPUT OF LOCAL ANCESTRY PIPELINE
        curr = pd.read_csv(tract_dir + f, delimiter='\t', names=["Chrom", "Genomic_Start", "Genomic_Stop", "Anc", "Start", "Stop"])
        curr["Info"] = ind_ID + "_" + hapl
        curr = curr.drop(["Genomic_Start", "Genomic_Stop"], axis=1)
        df = pd.merge(df, curr, how='outer')
    return(df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tract_dir', help='directory with local ancestry tracts')
    parser.add_argument('--out', help='output file')
    args = parser.parse_args()

    tracts_df = parse_tract_files(args.tract_dir)
    tracts_df = tracts_df[tracts_df.Anc == 'YRI']
    tracts_df = tracts_df.drop(["Anc"], axis=1)

    tracts_df.to_csv(args.out, sep='\t', index=False, header=False)

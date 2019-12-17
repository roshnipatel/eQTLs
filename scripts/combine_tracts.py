# Combines multiple local ancestry files (based on output from LocalAncestry
# pipeline) into a single BED-formatted file. Important: output is NOT sorted.

import argparse
import pandas as pd
import os
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--tract_dir', help='directory with local ancestry tracts')
parser.add_argument('--out', help='output file')
args = parser.parse_args()

tract_files = os.listdir(args.tract_dir)
tracts_df = pd.DataFrame(columns=["Chrom", "Start", "Stop", "Info"])
for f in tract_files:
    ind_ID = f.split('.')[2] # WILL NEED TO CHANGE THIS IF YOU UPDATE THE OUTPUT OF LOCAL ANCESTRY PIPELINE
    hapl = f.split('.')[3] # WILL NEED TO CHANGE THIS IF YOU UPDATE THE OUTPUT OF LOCAL ANCESTRY PIPELINE
    curr = pd.read_csv(args.tract_dir + f, delimiter='\t', names=["Chrom", "Genomic_Start", "Genomic_Stop", "Anc", "Start", "Stop"])
    curr = curr[curr.Anc == 'YRI']
    curr["Info"] = ind_ID + "_" + hapl
    curr = curr.drop(["Genomic_Start", "Genomic_Stop", "Anc"], axis=1)
    tracts_df = pd.merge(tracts_df, curr, how='outer')

tracts_df.to_csv(args.out, sep='\t', index=False, header=False)

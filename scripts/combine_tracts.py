# Combines multiple local ancestry files (based on output from LocalAncestry
# pipeline) into a single BED-formatted file. Important: output is NOT sorted.

import argparse
import pandas as pd
import os
import numpy as np

def parse_tract_files(tract_dir):
    all_bed_files = []
    for i in range(1, 23):
        chr_path = os.path.join(tract_dir, "chr{0}".format(str(i)))
        chr_files = [os.path.join(chr_path, f) for f in os.listdir(chr_path) if f[-3:] == "bed"]
        all_bed_files.extend(chr_files)
    df = pd.DataFrame()
    for f in all_bed_files:
        ind_id = f.split('/')[-1].split('.')[0] # WILL NEED TO CHANGE THIS IF YOU UPDATE THE OUTPUT OF LOCAL ANCESTRY PIPELINE
        hapl = f.split('/')[-1].split('.')[1] # WILL NEED TO CHANGE THIS IF YOU UPDATE THE OUTPUT OF LOCAL ANCESTRY PIPELINE
        curr = pd.read_csv(f, delimiter='\t', names=["chrom", "genomic_start", "genomic_stop", "anc", "start", "stop"])
        curr["info"] = curr.apply(lambda row: ind_id + "_" + hapl + "_" + row.anc, axis=1)
        curr = curr[["chrom", "start", "stop", "info"]]
        df = pd.concat([df, curr])
    return(df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tract_dir', help='directory with local ancestry tracts')
    parser.add_argument('--out', help='output file')
    args = parser.parse_args()

    tracts_df = parse_tract_files(args.tract_dir)
    tracts_df.to_csv(args.out, sep='\t', index=False, header=False)

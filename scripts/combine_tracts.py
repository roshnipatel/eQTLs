import argparse
import pandas as pd
import os
import numpy as np

def parse_tract_files(tract_dir):
    """LocalAncestry pipeline produces BED files of local ancestry tracts for
       each individual for each chromosome. This function loops over all files,
       combining all tracts into a single DataFrame. Function is hard-coded to 
       the exact output of the LocalAncestry pipeline, sorrynotsorry ¯\_(ツ)_/¯"""
    # Create list of all BED files in provided directory.
    all_bed_files = []
    for i in range(1, 23):
        chr_path = os.path.join(tract_dir, "chr{0}".format(str(i)))
        chr_files = [os.path.join(chr_path, f) for 
                     f in os.listdir(chr_path) if f[-3:] == "bed"]
        all_bed_files.extend(chr_files)

    # Parse each file and store to DataFrame. Again, this is HARD-CODED to the
    # current output of the LocalAncestry pipeline.
    df = pd.DataFrame()
    col_names = ["chrom", "genomic_start", "genomic_stop", "anc", "start", "stop"]
    for f in all_bed_files:
        ind_id = f.split('/')[-1].split('.')[0] 
        hapl = f.split('/')[-1].split('.')[1] 
        curr = pd.read_csv(f, delimiter='\t', names=col_names)
        metadata = ind_id + "_" + hapl + "_" + row.anc 
        curr["info"] = curr.apply(lambda row: metadata, axis=1)
        curr = curr[["chrom", "start", "stop", "info"]]
        df = pd.concat([df, curr])
    return(df)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--tract_dir')
    parser.add_argument('--out')
    args = parser.parse_args()

    tracts_df = parse_tract_files(args.tract_dir)
    tracts_df.to_csv(args.out, sep='\t', index=False, header=False)

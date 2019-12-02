#!/usr/bin/env python3
# Author: Francois Aguet
import pandas as pd
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description='Combine covariates into a single matrix')
parser.add_argument('expression_covariates', help='')
parser.add_argument('prefix', help='')
parser.add_argument('genotype_pcs', default=None, help='Genotype PCs')
parser.add_argument('samp_covariates', default=[], help='Additional covariates')
parser.add_argument('-o', '--output_dir', default='.', help='Output directory')
args = parser.parse_args()

print('Combining covariates ... ', end='', flush=True)
additional_df = pd.read_csv(args.samp_covariates, sep='\t', index_col=0, dtype=str)

# Pull out the categorical variables
plate = additional_df["BroadUWPlate"]
additional_df = additional_df.drop("BroadUWPlate", axis=1)

# Combine all the quantitative variables
expression_df = pd.read_csv(args.expression_covariates, sep='\t', index_col=0, dtype=str).transpose()
genotype_df = pd.read_csv(args.genotype_pcs, sep=' ', index_col=0, header=None, dtype=str)
ID_overlap = set(additional_df.index).intersection(set(expression_df.index)).intersection(set(genotype_df.index))
combined_df = pd.concat([additional_df.loc[ID_overlap], expression_df.loc[ID_overlap], genotype_df.loc[ID_overlap]], axis=1)

# identify and drop colinear covariates
C = combined_df.astype(np.float64)
Q,R = np.linalg.qr(C-np.mean(C, axis=0))
colinear_ix = np.abs(np.diag(R)) < np.finfo(np.float64).eps * C.shape[1]
if np.any(colinear_ix):
    print('Colinear covariates detected:')
    for i in C.columns[colinear_ix]:
        print("  * dropped '{}'".format(i))
    combined_df = combined_df.loc[~colinear_ix]

# Merge the categorical variables back in
combined_df = pd.concat([combined_df, plate], axis=1, join="inner").transpose()

combined_df.to_csv(args.prefix+'.combined_covariates.txt', sep='\t')
print('done.')

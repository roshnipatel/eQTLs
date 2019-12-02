import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sample_data', help="RNASeq sample metadata file")
parser.add_argument('output', help="file with TOR ID and NWDID of selected samples")
args = parser.parse_args()

samples = pd.read_csv(args.sample_data, delimiter='\t')

# Filter for PBMC RNASeq data
samples = samples[samples.Specimen == 'PBMC']

# Remove samples with replicates because we don't know how to deal with them right now
samples = samples.groupby(['Exam', 'NWDID']).filter(lambda x: len(x) == 1)

# Select data from visit 5 when possible
samples_5 = samples[samples.Exam == 5]

# For individuals without data from visit, select data from visit 1
ID_dif = set(samples.NWDID).difference(samples_5.NWDID)
where_dif = samples.NWDID.isin(ID_dif)
samples_1 = samples[where_dif]

# Merge data from visit 1 and 5
samples = samples_1.append(samples_5, ignore_index=True)

# Write TOR ID and NWDID to file
samples[['TOR_ID', 'NWDID']].to_csv(args.output, sep='\t', index=False)

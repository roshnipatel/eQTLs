import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--samples", help="file mapping RNASeq samples to individual IDs")
parser.add_argument("--indiv_metadata", help="file with individual metadata info")
parser.add_argument("--sample_metadata", help="file with sample metadata info")
parser.add_argument("--out", help="output file")
args = parser.parse_args()

samples = pd.read_csv(args.samples, delimiter='\t')
indiv_metadata = pd.read_csv(args.indiv_metadata)
sample_metadata = pd.read_csv(args.sample_metadata, delimiter='\t')

samples = pd.merge(samples, sample_metadata[["TOR_ID", "BroadUWPlate", "Exam", "ageatExam"]])
samples = pd.merge(samples, indiv_metadata[["TOR_ID_rnaseq", "gender1"]], left_on="TOR_ID", right_on="TOR_ID_rnaseq")
samples = samples.drop("TOR_ID_rnaseq", axis=1).drop("TOR_ID", axis=1).drop_duplicates()

samples.to_csv(args.out, sep='\t', index=False)

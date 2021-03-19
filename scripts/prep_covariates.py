import pandas as pd
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from call_eqtls import parse_genotypes

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", help="file mapping RNASeq samples to individual IDs", nargs='+')
    parser.add_argument("--sample_metadata", help="file with sample metadata info")
    parser.add_argument("--global_anc")
    parser.add_argument("--pca")
    parser.add_argument("--out")
    args = parser.parse_args()
    
    # Read in all samples
    samples = pd.DataFrame()
    for filepath in args.samples:
        tmp = pd.read_csv(filepath, delimiter='\t')
        samples = pd.concat([samples, tmp])
    
    # Merge sex, ancestry, and batch covariates
    sample_metadata = pd.read_csv(args.sample_metadata, delimiter='\t')
    ganc = pd.read_csv(args.global_anc, sep='\t')
    ganc = ganc[["ID", "YRI_admixture"]].rename({"ID": "nwd_id", 
                                                "YRI_admixture": "global_ancestry"}, axis=1)
    samples = pd.merge(samples, sample_metadata[["tor_id", "nwd_id",
                                                 "seq_center", "exam", 
                                                 "age", "sex", "race"]])
    samples = pd.merge(samples, ganc, how='left')
    samples["global_ancestry"].fillna(0, inplace=True)
    
    # Convert categorical variables to numerical quantities
    samples["seq_center"] = samples.apply(lambda row: 1 if row.seq_center == "UW" else 0, axis=1)
    samples["exam"] = samples.apply(lambda row: 1 if row.exam == "1" else 0, axis=1)
    samples["race"] = samples.apply(lambda row: 1 if "Black" in row.race else 0, axis=1)
    samples["sex"] = samples.apply(lambda row: 1 if row.sex == "Female" else 0, axis=1)
    
    # Merge all the covariates together
    for filepath in args.pca:
        PC_df = pd.read_csv(filepath, index_col=0, sep='\t')
        samples = pd.merge(samples, PC_df, left_on="nwd_id", right_index=True)
    
    samples.to_csv(args.out, sep='\t', index=False)

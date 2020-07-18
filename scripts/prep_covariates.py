import pandas as pd
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from call_eqtls import parse_genotypes

parser = argparse.ArgumentParser()
parser.add_argument("--samples", help="file mapping RNASeq samples to individual IDs")
parser.add_argument("--indiv_metadata", help="file with individual metadata info")
parser.add_argument("--sample_metadata", help="file with sample metadata info")
parser.add_argument("--vcf", help="genotype VCF")
parser.add_argument("--exp", help="normalized expression file")
parser.add_argument("--out", help="output file")
args = parser.parse_args()

# Specify number of PCs
num_exp_PC = 10
num_geno_PC = 15

def pca(n, X):
    scaled = StandardScaler().fit_transform(X)
    PCs = PCA(n_components=n).fit_transform(scaled)
    PC_df = pd.DataFrame(PCs, columns=["PC" + str(i) for i in range(n)])
    return(PC_df)

# Prepare age, sex, and batch covariates
samples = pd.read_csv(args.samples, delimiter='\t')
indiv_metadata = pd.read_csv(args.indiv_metadata)
sample_metadata = pd.read_csv(args.sample_metadata, delimiter='\t')
samples = pd.merge(samples, sample_metadata[["TOR_ID", "BroadUWPlate", "Exam", "ageatExam"]])
samples["Plate_Prefix"] = samples.apply(lambda row: 1 if row.BroadUWPlate[:2] == "SK" else 0, axis=1)
samples = pd.merge(samples, indiv_metadata[["TOR_ID_rnaseq", "gender1"]], left_on="TOR_ID", right_on="TOR_ID_rnaseq")
samples = samples.drop("TOR_ID_rnaseq", axis=1).drop("TOR_ID", axis=1).drop_duplicates()

# Generate expression PCs
expression = pd.read_csv(args.exp, sep='\t')
expression = expression.set_index("gene_id")
expression = expression.drop(["#chr", "start", "end"], axis=1)
expression = expression.dropna() # Remove genes where 1+ individuals have NaN values
expression = expression.T # Convert to individual x gene matrix in order to generate PCs
exp_PC = pca(num_exp_PC, expression) # Generate PCs on entire dataset
exp_PC["NWDID"] = list(expression.index)

# Generate genotype PCs
genotypes = parse_genotypes(args.vcf)
genotypes = genotypes.T # Convert to individual x genotype matrix in order to generate PCs
geno_PC = pca(num_geno_PC, genotypes)
geno_PC["NWDID"] = list(genotypes.index)

# Merge all the covariates together
PC_df = pd.merge(exp_PC, geno_PC, on="NWDID", suffixes=("_exp", "_geno"))
samples = pd.merge(samples, PC_df, on="NWDID")

samples.to_csv(args.out, sep='\t', index=False)

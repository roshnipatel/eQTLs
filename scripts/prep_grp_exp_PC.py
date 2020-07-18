import pandas as pd
import argparse
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from call_eqtls import parse_genotypes
import statsmodels.api as sm

parser = argparse.ArgumentParser()
parser.add_argument("--afr", help="file mapping RNASeq samples to individual IDs for AA")
parser.add_argument("--eur", help="file mapping RNASeq samples to individual IDs for EA")
parser.add_argument("--exp", help="normalized expression file", nargs='+')
parser.add_argument("--out", help="output file")
args = parser.parse_args()

# Specify number of PCs
num_exp_PC = 10

def pca(n, X):
    scaled = StandardScaler().fit_transform(X)
    PCs = PCA(n_components=n).fit_transform(scaled)
    PC_df = pd.DataFrame(PCs, columns=["PC" + str(i) for i in range(n)])
    PC_df.index = list(X.index)
    return(PC_df)

afr_samples = pd.read_csv(args.afr, delimiter='\t', index_col="NWDID")
eur_samples = pd.read_csv(args.eur, delimiter='\t', index_col="NWDID")

# Generate expression PCs
expression = pd.DataFrame()
expression_na_filled = pd.DataFrame()
join_type = "outer" # Need to have outer join for first merge with empty dataframe; second can be inner
for file_path in args.exp: 
    tmp = pd.read_csv(file_path, sep='\t')
    tmp = tmp.set_index("gene_id")
    tmp = tmp.drop(["#chr", "start", "end"], axis=1)
    tmp = tmp.T # Convert to individual x gene matrix in order to generate PCs
    expression = pd.concat([tmp, expression], join=join_type)
    tmp = tmp.fillna(tmp.median()) # Impute missing values with median - this is only used for calculating PCs
    expression_na_filled = pd.concat([tmp, expression_na_filled], join=join_type)
    join_type = "inner"

exp_PC = pca(num_exp_PC, expression_na_filled) # Generate PCs on entire dataset
exp_PC.to_csv(args.out + ".10_PC.txt", sep='\t')

afr_ind = list(afr_samples.index)
eur_ind = list(eur_samples.index)

def add_const(row):
    if row.name in afr_ind:
        return(pd.Series([1, 1, 0]))
    elif row.name in eur_ind:
        return(pd.Series([1, 0, 1]))
    else:
        return(pd.Series([1, 0, 0]))

# Add constants for regression
const = exp_PC.apply(add_const, axis=1) 
const.columns = ["global", "afr", "eur"]
cov = exp_PC.join(const)

def regress_covariates(row):
    model = sm.OLS(row, cov, missing='drop') # Drop missing values
    results = model.fit()
    resid_exp = results.resid
    return(resid_exp)

residual = expression.T.apply(regress_covariates, axis=1)
residual[afr_ind].to_csv(args.out + ".Afr.bed.gz", sep='\t')
residual[eur_ind].to_csv(args.out + ".Eur.bed.gz", sep='\t')

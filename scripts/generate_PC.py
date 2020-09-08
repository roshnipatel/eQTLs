import pandas as pd
import numpy as np
import argparse
from call_eqtls import parse_genotypes
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

parser = argparse.ArgumentParser()
parser.add_argument("--data", help="data file", nargs='+')
parser.add_argument("--n", help="number of PCs to report (if already known from e.g. scree plot)", default=None)
parser.add_argument("--filetype", help="either contains 'vcf' or 'expression'; parsed by Snakemake from specified file input")
parser.add_argument("--out", help="prefix of output files")
args = parser.parse_args()

# Generate PCs, imputing missing values with ancestry-specific median
data_df = pd.DataFrame()
join_type = "outer" # Need to have outer join for first merge with empty dataframe; second can be inner
for file_path in args.data: 
    if "vcf" in args.filetype:
        filetype = "genotype"
        tmp = parse_genotypes(file_path)
    elif "expression" in args.filetype:
        filetype = "expression"
        tmp = pd.read_csv(file_path, sep='\t')
        tmp = tmp.set_index("gene_id")
        tmp = tmp.drop(["#chr", "start", "end"], axis=1)
    tmp = tmp.T # Convert to individual x variable matrix in order to generate PCs
    tmp = tmp.fillna(tmp.median()) 
    data_df = pd.concat([tmp, data_df], join=join_type)
    join_type = "inner"

scaled_df = StandardScaler().fit_transform(data_df)
PCA_obj = PCA().fit(scaled_df)

# Write eigenvalues to file
eigenvals = PCA_obj.explained_variance_
np.savetxt(args.out + ".eigenvalues.txt", eigenvals, delimiter='\n')

# Write components to file
if args.n is not None:
    num_PC = int(args.n)
    PCA_obj = PCA_obj.set_params(n_components = num_PC)
    components = PCA_obj.fit_transform(scaled_df)
    PC_df = pd.DataFrame(components, columns=[filetype + "_PC" + str(i) for i in range(num_PC)])
    PC_df.index = list(data_df.index)
    PC_df.to_csv(args.out + ".principal_components.txt", sep='\t')
else:
    components = PCA_obj.transform(scaled_df)
    PC_df = pd.DataFrame(components, columns=[filetype + "_PC" + str(i) for i in range(len(eigenvals))])
    PC_df.index = list(data_df.index)
    PC_df.to_csv(args.out + ".principal_components.txt", sep='\t')

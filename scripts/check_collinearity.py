import numpy as np
import pandas as pd

def check_eigenvals(group):
    pivoted_group = group[["NWDID", "Genotype_Afr", "Genotype_Eur", "Race_AA", "Gene"]].pivot(index="NWDID", columns="Gene", values=["Genotype_Afr", "Genotype_Eur", "Race_AA"])
    pivoted_group["Race_EA"] = 1 - pivoted_group["Race_AA"]
    design_matrix = np.array(pivoted_group.reset_index().iloc[:,1:])
    gram_matrix = np.matmul(design_matrix.T, design_matrix)
    eigvals = np.linalg.eigvals(gram_matrix)
    all_positive_eigvals = all(v > 0.001 for v in eigvals)
    return(all_positive_eigvals)

merged = pd.read_csv("data/MLE/merged_data.txt", sep='\t')
eigval_sign = merged.groupby("Gene").apply(check_eigenvals)
all_true = all(v for v in eigval_sign)
print(all_true)

import pandas as pd
import argparse
import numpy as np
import cvxpy as cp
import statsmodels.api as sm

def document_params(delta_file, betas_file, delta, betas):
    with open(betas_file, 'w') as f:
        betas.to_csv(f, sep='\t')
    with open(delta_file, 'a') as f:
        f.write(str(delta))
        f.write('\n')

def optimize_delta(df):
    b = np.array(df["Expression"] - df["Genotype_Afr"] * df["Curr_Effect_Afr"] - df["Genotype_Eur"] * df["Curr_Effect_Eur"] - df["Curr_Int_Afr"] * df["Race_AA"] - df["Curr_Int_Eur"] * (1 - df["Race_AA"]))
    A = np.array((df["Curr_Effect_Afr"] - df["Curr_Effect_Eur"]) * df["Genotype_Eur"] * df["Race_AA"])
    P = np.matrix(np.dot(A, A))
    q = np.matrix(np.dot(A, -b))
    G = np.array([1, -1])
    h = np.array([1, 0])
    delta = cp.Variable(1)
    prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(delta, P) + q.T @ delta), [G @ delta <= h])
    prob.solve()
    result = delta.value[0]
    return(result)

def optimize_betas(group):
    X_Afr = group["Genotype_Afr"] + group["Curr_Delta"] * group["Genotype_Eur"] * group["Race_AA"]
    X_Eur = group["Genotype_Eur"] - group["Curr_Delta"] * group["Genotype_Eur"] * group["Race_AA"]
    Int_Afr = group["Race_AA"]
    Int_Eur = 1 - group["Race_AA"]
    Y = group["Expression"]
    # df = pd.DataFrame({"Y": Y, "X_Afr": X_Afr, "X_Eur": X_Eur})
    # result = sm.OLS(df["Y"], df[["X_Afr", "X_Eur"]]).fit().params
    df = pd.DataFrame({"Y": Y, "X_Afr": X_Afr, "X_Eur": X_Eur, "Int_Afr": Int_Afr, "Int_Eur": Int_Eur})
    result = sm.OLS(df["Y"], df[["X_Afr", "X_Eur", "Int_Afr", "Int_Eur"]]).fit().params
    return(result)

def update_params(df, betas=None, delta=None):
    if delta is not None:
        df["Curr_Delta"] = delta
    elif betas is not None:
        # betas = betas.rename(columns={"X_Afr": "Curr_Effect_Afr", "X_Eur": "Curr_Effect_Eur"})
        # df = df.drop(columns=["Curr_Effect_Afr", "Curr_Effect_Eur"])
        betas = betas.rename(columns={"X_Afr": "Curr_Effect_Afr", "X_Eur": "Curr_Effect_Eur", "Int_Afr": "Curr_Int_Afr", "Int_Eur": "Curr_Int_Eur"})
        df = df.drop(columns=["Curr_Effect_Afr", "Curr_Effect_Eur", "Curr_Int_Afr", "Curr_Int_Eur"])
        df = pd.merge(df, betas, left_on="Gene", right_index=True)
    return(df)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--max_iter')
    parser.add_argument('--delta_out')
    parser.add_argument('--betas_out')
    args = parser.parse_args()
    
    merged_data = pd.read_csv(args.merged, sep='\t')
    
    merged_data[["Curr_Effect_Afr", "Curr_Effect_Eur"]] = merged_data[["effect_Afr", "effect_Eur"]]
    merged_data[["Curr_Int_Afr", "Curr_Int_Eur"]] = merged_data[["intercept_Afr", "intercept_Eur"]]
    prev_delta = -1 # Ensures we do 1+ iterations
    for i in range(int(args.max_iter)):
        curr_delta = optimize_delta(merged_data)
        if (abs(curr_delta - prev_delta) < .0001):
            break
        prev_delta = curr_delta
        merged_data = update_params(merged_data, delta=curr_delta)
        curr_betas = merged_data.groupby("Gene").apply(optimize_betas)
        merged_data = update_params(merged_data, betas=curr_betas)
        document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)

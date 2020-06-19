import pandas as pd
import argparse
import numpy as np
import cvxpy as cp
import statsmodels.api as sm

def document_params(delta_file, betas_file, delta, betas):
    """Write current iteration of parameters to file."""
    with open(betas_file, 'w') as f:
        betas.to_csv(f, sep='\t')
    with open(delta_file, 'a') as f:
        f.write(str(delta))
        f.write('\n')

def optimize_delta(df):
    """Optimize delta in likelihood model using quadratic programming solver."""
    # # Optimize delta based on initial likelihood model (estimating African and European effects separately)
    # b = np.array(df["Expression"] - df["Genotype_Afr"] * df["Curr_Effect_Afr"] - df["Genotype_Eur"] * df["Curr_Effect_Eur"] - df["Curr_Int_Afr"] * df["Race_AA"] - df["Curr_Int_Eur"] * (1 - df["Race_AA"]))
    # A = np.array((df["Curr_Effect_Afr"] - df["Curr_Effect_Eur"]) * df["Genotype_Eur"] * df["Race_AA"])
    # h = np.array([1, 0])

    # Optimize delta based on new likelihood model (estimating mean effect and half the difference between African and European effects)
    b = np.array(df["Expression"] - (df["Genotype_Afr"] + df["Genotype_Eur"]) * df["Curr_Effect_mean"] - df["Genotype_Afr"] * df["Curr_Effect_Afr"] + df["Genotype_Eur"] * df["Curr_Effect_Afr"] - df["Curr_Int_Afr"] * df["Race_AA"] - df["Curr_Int_Eur"] * (1 - df["Race_AA"]))
    A = np.array(df["Genotype_Eur"] * df["Curr_Effect_Afr"] * df["Race_AA"])
    h = np.array([2, 0])
    
    # Set up quadratic programming problem and solve
    P = np.matrix(np.dot(A, A))
    q = np.matrix(np.dot(A, -b))
    G = np.array([1, -1])
    delta = cp.Variable(1)
    prob = cp.Problem(cp.Minimize((1/2)*cp.quad_form(delta, P) + q.T @ delta), [G @ delta <= h])
    prob.solve()
    result = delta.value[0]
    return(result)

def optimize_betas(group):
    """Optimize betas/effects in likelihood model using ordinary least squares linear regression."""
    # # Optimize betas (effects) based on initial likelihood model (estimating African and European effects separately)
    # X_Afr = group["Genotype_Afr"] + group["Curr_Delta"] * group["Genotype_Eur"] * group["Race_AA"]
    # X_Eur = group["Genotype_Eur"] - group["Curr_Delta"] * group["Genotype_Eur"] * group["Race_AA"]
    # Int_Afr = group["Race_AA"]
    # Int_Eur = 1 - group["Race_AA"]
    # Y = group["Expression"]
    # df = pd.DataFrame({"Y": Y, "X_Afr": X_Afr, "X_Eur": X_Eur, "Int_Afr": Int_Afr, "Int_Eur": Int_Eur})
    # result = sm.OLS(df["Y"], df[["X_Afr", "X_Eur", "Int_Afr", "Int_Eur"]]).fit().params
    
    # Optimize betas (effects) based on new likelihood model (estimating mean effect and half the difference between African and European effects)
    X_mean = group["Genotype_Afr"] + group["Genotype_Eur"]
    X_Afr = group["Genotype_Afr"] - group["Genotype_Eur"] + group["Curr_Delta"] * group["Genotype_Eur"] * group["Race_AA"]
    Int_Afr = group["Race_AA"]
    Int_Eur = 1 - group["Race_AA"]
    Y = group["Expression"]
    df = pd.DataFrame({"Y": Y, "X_Afr": X_Afr, "X_mean": X_mean, "Int_Afr": Int_Afr, "Int_Eur": Int_Eur})
    result = sm.OLS(df["Y"], df[["X_Afr", "X_mean", "Int_Afr", "Int_Eur"]]).fit().params

    return(result)

def update_params(df, betas=None, delta=None):
    """Update current parameters in provided dataframe."""
    if delta is not None:
        df["Curr_Delta"] = delta
    elif betas is not None:
        # # Update parameter names based on initial likelihood model (estimating African and European effects separately)
        # betas = betas.rename(columns={"X_Afr": "Curr_Effect_Afr", "X_Eur": "Curr_Effect_Eur", "Int_Afr": "Curr_Int_Afr", "Int_Eur": "Curr_Int_Eur"})
        # df = df.drop(columns=["Curr_Effect_Afr", "Curr_Effect_Eur", "Curr_Int_Afr", "Curr_Int_Eur"])

        # Update parameter names based on new likelihood model (estimating mean effect and half the difference between African and European effects) 
        betas = betas.rename(columns={"X_Afr": "Curr_Effect_Afr", "X_mean": "Curr_Effect_mean", "Int_Afr": "Curr_Int_Afr", "Int_Eur": "Curr_Int_Eur"})
        df = df.drop(columns=["Curr_Effect_Afr", "Curr_Effect_mean", "Curr_Int_Afr", "Curr_Int_Eur"])

        df = pd.merge(df, betas, left_on="Gene", right_index=True)
    return(df)

def simulate_data(delta, afr_mean, eur_mean, afr_effect, eur_effect, afr_maf, eur_maf):
    """Simulate some mock data for troubleshooting optimization code."""
    def simulate_pheno(row):
        mean_effect = (afr_effect + eur_effect) / 2
        half_dif = (afr_effect - eur_effect) / 2
        genetic_effects = (row.Genotype_Eur + row.Genotype_Afr) * mean_effect + half_dif * row.Genotype_Afr - half_dif * row.Genotype_Eur
        dif_pop_mean = row.Race_AA * afr_mean + (1 - row.Race_AA) * eur_mean
        interaction = delta * half_dif * row.Genotype_Eur * row.Race_AA
        noise = np.random.normal()
        return(genetic_effects + dif_pop_mean + interaction + noise)
    n_hom = 100
    n_het = 50
    n_Eur = 380
    AA_hom_geno_Afr = np.random.binomial(2, afr_maf, n_hom)
    AA_hom_geno_Eur = np.zeros(n_hom)
    AA_het_geno_Afr = np.random.binomial(1, afr_maf, n_het)
    AA_het_geno_Eur = np.random.binomial(1, eur_maf, n_het)
    EA_geno_Eur = np.random.binomial(2, eur_maf, n_Eur)
    EA_geno_Afr = np.zeros(n_Eur)
    geno_Eur = np.hstack([AA_hom_geno_Eur, AA_het_geno_Eur, EA_geno_Eur])
    geno_Afr = np.hstack([AA_hom_geno_Afr, AA_het_geno_Afr, EA_geno_Afr])
    race_AA = np.hstack([np.ones(n_hom + n_het), np.zeros(n_Eur)])
    merged_df = pd.DataFrame({"Genotype_Eur": geno_Eur, "Genotype_Afr": geno_Afr, "Race_AA": race_AA})
    merged_df["Expression"] = merged_df.apply(simulate_pheno, axis=1)
    merged_df["Gene"] = "Test"
    return(merged_df)

def neg_control(group):
    """Filter for individuals in African and European estimation sets. Code Validation Eur as AA in order to 
       perform a negative control. Optimizing model on this recoded dataset should result in delta = 0."""
    gene = group.name
    Afr_idv = []
    with open("data/fastqtl_sample_input/estimation/Afr/" + gene + ".txt", 'r') as f:
        for idv in f:
            Afr_idv.append(idv.strip())
    Eur_idv = []
    with open("data/fastqtl_sample_input/estimation/Eur/" + gene + ".txt", 'r') as f:
        for idv in f:
            Eur_idv.append(idv.strip())
    val_idv = []
    with open("data/fastqtl_sample_input/validation/Eur/" + gene + ".txt", 'r') as f:
        for idv in f:
            val_idv.append(idv.strip())
    group_subset = group[(group.NWDID.isin(Afr_idv)) | (group.NWDID.isin(Eur_idv)) | (group.NWDID.isin(val_idv))]
    group_subset["Race_AA"] = group_subset.apply(lambda row: 1 if row.NWDID in val_idv else row.Race_AA, axis=1, result_type='expand')
    return(group_subset)

def drop_asc(group):
    """Remove European ascertainment individuals from dataset for optimizing likelihood model. These 
       individuals were used to select the genes/variants we are not optimizing the model on, so inclusion
       of these individuals will result in inflated European effect sizes."""
    gene = group.name
    asc_idv = []
    with open("data/fastqtl_sample_input/ascertainment/Eur/" + gene + ".txt", 'r') as f:
        for idv in f:
            asc_idv.append(idv.strip())
    return(group[-group.NWDID.isin(asc_idv)])

def bootstrap_data(all_df):
    """Generate one random sample (with replacement) from the input dataframe. Size of sample is 
       specified by number of genes in input dataframe."""
    def decrement(d): # No return statement because dictionaries are modified in-place
        keys_to_drop = []
        for key, val in d.items():
            val = val - 1
            if val == 0:
                keys_to_drop.append(key)
            else:
                d[key] = val
        for key in keys_to_drop:
            del d[key]

    gene_list = all_df["Gene"].unique()
    n_genes = len(gene_list)
    bootstrap_list = np.random.choice(gene_list, n_genes, replace=True)
    bootstrap_dict = {}
    for gene in bootstrap_list:
        bootstrap_dict[gene] = bootstrap_dict.get(gene, 0) + 1

    bootstrap_df = pd.DataFrame()
    while bootstrap_dict:
        tmp_df = all_df[all_df["Gene"].isin(bootstrap_dict.keys())]
        bootstrap_df = pd.concat([bootstrap_df, tmp_df])
        decrement(bootstrap_dict)

    return(bootstrap_df)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--merged')
    parser.add_argument('--max_iter')
    parser.add_argument('--delta_out')
    parser.add_argument('--betas_out')
    parser.add_argument('--option')
    parser.add_argument('--rand_start', action='store_true')
    parser.add_argument('--simulate_data', action='store_true')
    parser.add_argument('--delta', default=None)
    parser.add_argument('--bootstrap', action='store_true')
    parser.add_argument('--sim_afr_mean')
    parser.add_argument('--sim_eur_mean')
    parser.add_argument('--sim_afr_effect')
    parser.add_argument('--sim_eur_effect')
    parser.add_argument('--sim_afr_maf')
    parser.add_argument('--sim_eur_maf')
    parser.add_argument('--sim_delta')
    args = parser.parse_args()
    
    if not args.simulate_data: # Read in genotypes/phenotypes from file
        merged_data = pd.read_csv(args.merged, sep='\t')
        
        # # Initialize params based on initial likelihood model (estimating African and European effects separately)
        # merged_data[["Curr_Effect_Afr", "Curr_Effect_Eur"]] = merged_data[["effect_Afr", "effect_Eur"]]
        # merged_data[["Curr_Int_Afr", "Curr_Int_Eur"]] = merged_data[["intercept_Afr", "intercept_Eur"]]
        
        # Initialize params based on new likelihood model (estimating mean effect and half the difference between African and European effects)
        merged_data["Curr_Effect_mean"] = (merged_data["effect_Afr"] + merged_data["effect_Eur"]) / 2
        merged_data["Curr_Effect_Afr"] = (merged_data["effect_Afr"] - merged_data["effect_Eur"]) / 2
        merged_data[["Curr_Int_Afr", "Curr_Int_Eur"]] = merged_data[["intercept_Afr", "intercept_Eur"]]

        # Deal with command line arguments/options
        if args.option == "neg_control":
            merged_data = merged_data.groupby("Gene").apply(neg_control).reset_index(drop=True)
        if args.option == "drop_asc":
            merged_data = merged_data.groupby("Gene").apply(drop_asc).reset_index(drop=True)
        if args.bootstrap:
            merged_data = bootstrap_data(merged_data)
    else: # Simulate genotypes/phenotypes
        merged_data = simulate_data(float(args.sim_delta), float(args.sim_afr_mean), float(args.sim_eur_mean), float(args.sim_afr_effect), float(args.sim_eur_effect), float(args.sim_afr_maf), float(args.sim_eur_maf))
        # # Initialize params based on initial likelihood model (estimating African and European effects separately)
        # merged_data["Curr_Effect_Afr"] = None
        # merged_data["Curr_Effect_Eur"] = None
        
        # Initialize params based on new likelihood model (estimating mean effect and half the difference between African and European effects)
        merged_data["Curr_Effect_Afr"] = None
        merged_data["Curr_Effect_mean"] = None
        merged_data["Curr_Int_Afr"] = None
        merged_data["Curr_Int_Eur"] = None

    if args.delta is not None: # Optimize betas for a user-specified delta + write to file
        curr_delta = args.delta
        merged_data = update_params(merged_data, delta=curr_delta)
        curr_betas = merged_data.groupby("Gene").apply(optimize_betas)
        merged_data = update_params(merged_data, betas=curr_betas)
        document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)
    else: # Iteratively optimize delta + betas, using linear regression estimates of betas as a jumping-off point
        prev_delta = -1 # Ensures we do 1+ iterations
        for i in range(int(args.max_iter)):
            if i == 0 and (args.rand_start or args.simulate_data):
                curr_delta = np.random.uniform(0, 1)
            else:
                curr_delta = optimize_delta(merged_data)
            if (abs(curr_delta - prev_delta) < .0001):
                break
            prev_delta = curr_delta
            merged_data = update_params(merged_data, delta=curr_delta)
            curr_betas = merged_data.groupby("Gene").apply(optimize_betas)
            merged_data = update_params(merged_data, betas=curr_betas)
            document_params(args.delta_out, args.betas_out, curr_delta, curr_betas)

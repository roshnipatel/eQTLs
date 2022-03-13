import stdpopsim
import pandas as pd
import tskit as tsk
import numpy as np
from numba import njit, prange
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--ganc')
parser.add_argument('--gene')
parser.add_argument('--delta')
parser.add_argument('--asc_pop')
parser.add_argument('--same_beta', action="store_true")
parser.add_argument('--r2')
parser.add_argument('--out', help="output")
args = parser.parse_args()

# Simulate number of local African haplotypes from global African ancestry frac
ganc = np.loadtxt(args.ganc)
lanc = np.random.binomial(2, ganc)

N_AA = lanc.shape[0] # Number of Afr-Am ind to simulate
N_EA = 500 # Number of Eur-Am ind to simulate
N_AFR = lanc.sum() # Number of African genomes to simulate
N_EUR = (N_AA * 2 - N_AFR) + (N_EA * 2) # Number of European genomes to simulate

# Specify number of ascertainment individuals depending on population
if args.asc_pop == "Eur":
    N_ASC = 250 
else:
    N_ASC = round((lanc == 2).sum() * .66)

# Simulate African + European genomes
species = stdpopsim.get_species("HomSap")
contig = species.get_contig("chr22", length_multiplier=0.01)
model = species.get_demographic_model("OutOfAfrica_2T12")
samples = model.get_samples(N_AFR, N_EUR)
engine = stdpopsim.get_engine('msprime')
ts = engine.simulate(model, contig, samples)
gts = ts.genotype_matrix()

# Apply MAF threshold for African and European ancestries separately
gts_afr = gts[:, 0:N_AFR]
gts_eur = gts[:, N_AFR:]
afr_freq = gts_afr.mean(axis=1)
eur_freq = gts_eur.mean(axis=1)
gts = gts[np.logical_and(np.logical_and(eur_freq < .95, eur_freq > .05),
                         np.logical_and(afr_freq < .95, afr_freq > .05)), :]
gts_afr = gts[:, 0:N_AFR]
gts_eur = gts[:, N_AFR:]

# Simulate causal effects
if args.same_beta:
    causal_eff = np.random.normal(0, 1)
    causal_eff = np.array([[causal_eff, causal_eff]])
else:
    var_cov = np.array([[1, 0.85], [0.85, 1]])
    mean = np.array([0, 0])
    causal_eff = np.random.multivariate_normal(mean, var_cov, 1)

@njit(parallel=True)
def collapse_anc_hom(mat):
    summed = np.zeros((mat.shape[0], mat.shape[1] // 2))
    for i in prange(mat.shape[1]):
        summed[:, i//2] += mat[:, i]
    return(summed)

# Simulate causal and tag SNP by jointly sampling from the set of all
# pairs of SNPs with r2 greater than the specified threshold. If 
# simulating under a regime of identical causal effects, set tag SNP 
# equal to causal SNP with probability 0.25.
if args.asc_pop == "Eur":
    geno_asc = collapse_anc_hom(gts_eur[:, :N_ASC * 2])
else:
    geno_asc = collapse_anc_hom(gts_afr[:, :N_ASC * 2])
corr = np.corrcoef(geno_asc)
r2_threshold = float(args.r2)
linked = np.argwhere(corr > r2_threshold)
idx = np.random.randint(0, linked.shape[0])
causal_idx = linked[idx, 0]
tag_idx = linked[idx, 1]
if args.same_beta:
    if np.random.rand() < 0.5:
        tag_idx = causal_idx

# Collapse causal and tag genotypes per-individual (only relevant for
# ancestry-homozygous individuals)
if args.asc_pop == "Eur":
    n_eur_sampled, n_afr_sampled = N_ASC * 2, 0
else:
    n_afr_sampled, n_eur_sampled = N_ASC * 2, 0

n_aa_0 = (lanc == 0).sum()
samp = n_aa_0 * 2
geno_aa_0 = collapse_anc_hom(gts_eur[[causal_idx, tag_idx], n_eur_sampled:n_eur_sampled + samp])
n_eur_sampled += samp

n_aa_1 = (lanc == 1).sum()
samp = n_aa_1
geno_aa_1_afr = gts_afr[[causal_idx, tag_idx], n_afr_sampled:n_afr_sampled + samp]
geno_aa_1_eur = gts_eur[[causal_idx, tag_idx], n_eur_sampled:n_eur_sampled + samp]
n_afr_sampled += samp
n_eur_sampled += samp

n_aa_2 = (lanc == 2).sum()
if args.asc_pop == "Eur":
    samp = n_aa_2 * 2
else:
    samp = (n_aa_2 - N_ASC) * 2
geno_aa_2 = collapse_anc_hom(gts_afr[[causal_idx, tag_idx], n_afr_sampled:n_afr_sampled + samp])
n_afr_sampled += samp

if args.asc_pop == "Eur":
    samp = (N_EA - N_ASC) * 2
else:
    samp = N_EA * 2
geno_ea = collapse_anc_hom(gts_eur[[causal_idx, tag_idx], n_eur_sampled:n_eur_sampled + samp])
n_eur_sampled += samp

# Reindex causal and tag SNP based on new geno matrices
causal_idx, tag_idx = 0, 1

# Simulate phenotypes from causal effects
delta = float(args.delta)
bA = causal_eff[0, 0]
bE = causal_eff[0, 1]
pheno_aa_0 = geno_aa_0[causal_idx, :] * (bA * delta + \
                                         bE * (1 - delta))
pheno_aa_1 = geno_aa_1_afr[causal_idx, :] * bA + \
             geno_aa_1_eur[causal_idx, :] * (bA * delta + \
                                             bE * (1 - delta))
pheno_aa_2 = geno_aa_2[causal_idx, :] * bA
pheno_ea = geno_ea[causal_idx, :] * bE

# Write phenos + anc-phased geno for tag SNP to file
df_aa_0 = pd.DataFrame({"expression": pheno_aa_0, "genotype_Eur": geno_aa_0[tag_idx, :],
                        "genotype_Afr": np.zeros(pheno_aa_0.shape[0]), "race": 1,
                        "gene": args.gene})

df_aa_1 = pd.DataFrame({"expression": pheno_aa_1, "genotype_Eur": geno_aa_1_eur[tag_idx, :],
                        "genotype_Afr": geno_aa_1_afr[tag_idx, :], "race": 1,
                        "gene": args.gene})

df_aa_2 = pd.DataFrame({"expression": pheno_aa_2, "genotype_Afr": geno_aa_2[tag_idx, :],
                        "genotype_Eur": np.zeros(pheno_aa_2.shape[0]), "race": 1,
                        "gene": args.gene})

df_ea = pd.DataFrame({"expression": pheno_ea, "genotype_Eur": geno_ea[tag_idx, :],
                      "genotype_Afr": np.zeros(pheno_ea.shape[0]), "race": 0,
                      "gene": args.gene})

comb_df = pd.concat([df_aa_0, df_aa_1, df_aa_2, df_ea], ignore_index=True)
comb_df.index.rename("nwd_id", inplace=True)
comb_df.to_csv(args.out, index=True, sep='\t')

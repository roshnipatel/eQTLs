import pandas as pd
import numpy as np

### Directories ###

# Directory in which to find protected data (anything linking NWDIDs to data)
PROTECTED_DATA_DIR = "protected_data/"

# Directory in which to find remaining data
DATA_DIR = "data/"

# Directory in which to store results (eQTL hits, plots, etc)
RES_DIR = "results/"



### Data ###

# Gencode annotation file (downloaded from gencodegenes.org/human/releases.html)
GENCODE_ANNO = "gencode.v30.annotation.gtf"

# Chromosome length file
CHR_LEN = "chrom_lengths.tsv"

# Genotype VCF
VCF = "genotype_freeze.6a.pass_only.phased.mesa_1331samples.maf01.biallelic.vcf.gz"

# Paths to RNASeq data
READS = "TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_reads.gct.gz"
TPM = "TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_tpm.gct.gz"

# Paths to metadata files
RNASEQ_METADATA = "metadata-TOPMed_MESA_RNAseq_2973samples_metadata.tsv"

# Paths to files from local and global ancestry determination
EXCLUSION_FILE = "individuals_to_exclude.txt" # File of individual IDs that failed QC
BED_DIR = "bed/" # Directory of local ancestry tracts



### Scripts ###

ANNO_SCRIPT = "scripts/collapse_annotation.py"
SAMPLE_SELECTION_SCRIPT = "scripts/select_samples.py"
NORM_SCRIPT = "scripts/normalize_expression.py"
COV_PREP_SCRIPT = "scripts/prep_covariates.py"
TRACT_SCRIPT = "scripts/combine_tracts.py"
GENE_BED_SCRIPT = "scripts/make_genes_bed.py"
PARTITION_SCRIPT = "scripts/partition_samples.py"
ESTIMATION_SCRIPT = "scripts/estimate_effect_sizes.py"
ID_HITS_SCRIPT = "scripts/identify_hits.py"
GENE_LIST_SCRIPT = "scripts/merge_gene_list.py"
CONCAT_QTL_SCRIPT = "scripts/concatenate_estimates.py"
PCA_SCRIPT = "scripts/generate_PC.py"
MERGE_SCRIPT = "scripts/merge_eQTL_data.py"
OPTIMIZE_SCRIPT = "scripts/iterative_parameter_optimization.py"
LIKELIHOOD_SCRIPT = "scripts/compute_likelihood.py"
EXP_CORRELATION_SCRIPT = "scripts/correlate_expression_covariates.py"
ANOVA_SCRIPT = "scripts/anova.py"
BOOTSTRAP_SCRIPT = "scripts/parse_bootstrap.py"
MASK_AA_GENO_SCRIPT = "scripts/mask_AA_genotypes.py"
MASK_EA_GENO_SCRIPT = "scripts/mask_EA_genotypes.py"
SIM_SCRIPT = "scripts/stdpopsim_simulations.py"



### Parameters ###

# Radius of cis window for testing SNPs
WINDOW = 100000 # 100 kb

# Number of iterations for likelihood EM
MAX_ITER = 300

# FDR threshold used
FDR = .01

# Chromosomes
CHROMS = [str(i) for i in range(1, 23)]

# Covariates
sites = "site1c_3.site1c_4.site1c_5.site1c_6.site1c_7.site1c_8"
FULL_COV_LIST = ["local_ancestry.global_ancestry.genotype_PC1.race_Afr.race_Eur.seq_center.exam.sex." + sites]
NO_BETA_COV_LIST = ["seq_center.exam." + sites, "seq_center.exam.sex." + sites,
                    "race_Afr.race_Eur.seq_center.exam.sex." + sites,
                    "global_ancestry.genotype_PC1.race_Afr.race_Eur.seq_center.exam.sex." + sites] + FULL_COV_LIST
GROUPS = ["test", "control"]



### Functions ###

def fetch_genes(f):
    genes, gene_idx = [], None
    for line in f.readlines():
        if not gene_idx:
            gene_idx = line.strip().split('\t').index("gene")
            continue
        curr_gene = line.strip().split('\t')[gene_idx]
        genes.append(curr_gene)
    return(genes)

import pandas as pd

### Directories ###

# Directory in which to find WGS and RNASeq data
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
EXCLUSION_FILE = "data/individuals_to_exclude.txt" # File of individual IDs that failed QC
BED_DIR = "data/bed/" # Directory of local ancestry tracts



### Scripts ###

ANNO_SCRIPT = "scripts/collapse_annotation.py"
SAMPLE_SELECTION_SCRIPT = "scripts/select_samples.py"
NORM_SCRIPT = "scripts/normalize_expression.py"
PEER_SCRIPT = "scripts/src/run_PEER.R"
COV_PREP_SCRIPT = "scripts/prep_covariates.py"
COV_COMB_SCRIPT = "scripts/src/combine_covariates.py"
TRACT_SCRIPT = "scripts/combine_tracts.py"
GENE_BED_SCRIPT = "scripts/make_genes_bed.py"
PARTITION_SCRIPT = "scripts/partition_samples.py"
EQTL_SCRIPT = "scripts/call_eqtls.py"
ID_HITS_SCRIPT = "scripts/identify_hits.py"
MERGE_GENES = "scripts/merge_gene_list.py"
MERGE_RES = "scripts/merge_eqtl_calls.py"
PCA_SCRIPT = "scripts/generate_PC.py"
SD_SCRIPT = "scripts/std_dev_comparison.R"
SIM_SCRIPT = "scripts/simulation.py"



### Misc. Utilities ###

# Docker environment containing FastQTL (downloaded from https://hub.docker.com/r/broadinstitute/gtex_eqtl/)
FASTQTL_DOCKER = "gtex_eqtl_V8.sif"



### Parameters ###

# Window size for eQTL calling
WINDOW = 100000 # 100 kb

# Permutation range for fastQTL permutation pass
PERM = "1000 1000"

# Number of iterations for likelihood EM
MAX_ITER = 300

# FDR threshold used
FDR = .05

# Chromosomes
CHROMS = [str(i) for i in range(1, 23)]

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
RNASEQ_METADATA = "MESA_TOPMed_RNASeqSamples_11022018.txt"
IND_METADATA = "MESA_sample_info.csv"



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
HITS_SCRIPT = "scripts/identify_hits.R"
MERGE_GENES = "scripts/merge_gene_list.py"
MERGE_RES = "scripts/merge_eqtl_calls.py"


### Misc. Utilities ###

# Docker environment containing FastQTL (downloaded from https://hub.docker.com/r/broadinstitute/gtex_eqtl/)
FASTQTL_DOCKER = "gtex_eqtl_V8.sif"



### Parameters ###

# Window size for eQTL calling
WINDOW = 100000 # 100 kb

# Permutation range for fastQTL permutation pass
PERM = "1000 1000"


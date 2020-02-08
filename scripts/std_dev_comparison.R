library(tidyverse)
library(matrixStats)

reads = read_tsv("data/TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_reads.gct.gz", skip=2)
metadata = read_csv("data/MESA_sample_info.csv")
samples = read_tsv("data/sample_participant_lookup.txt")
hits = read_tsv("results/v2/hits/hits_ascertainment_Eur.txt")
filt_genes = read_tsv("data/mesa.expression.bed.gz")

# Filter for AA and EA samples by race information in individual metadata
AA = samples %>% left_join(metadata %>% select(NWDID, race1c) %>% unique) %>% filter(race1c == 3)
EA = samples %>% left_join(metadata %>% select(NWDID, race1c) %>% unique) %>% filter(race1c == 1)
# EA = EA %>% sample_n(AA %>% nrow) # Optionally downsample EA dataset until it has the same sample size as AA
AAreads = reads %>% select(AA$TOR_ID)
EAreads = reads %>% select(EA$TOR_ID)

# Within each population, compute standard deviation in gene reads for all genes
AAsd = AAreads %>% transmute(AA_sd = rowSds(as.matrix(.[AAreads %>% colnames])))
EAsd = EAreads %>% transmute(EA_sd = rowSds(as.matrix(.[EAreads %>% colnames])))

# Normalize expression SD measurement relative to AA expression SD
genes = bind_cols(reads %>% select(Name, Description), AAsd, EAsd)
genes = genes %>% mutate(AA_norm = 1, EA_norm = EA_sd/AA_sd)
# genes %>% write.table("data/expression_std_dev.txt", row.names=FALSE, sep='\t')

# Filter for genes that passed pipeline filters (e.g. min # reads in min # ind) and significant eGenes
genes_filt = genes %>% inner_join(filt_genes %>% transmute(Name = gene_id))
genes_hits = genes_filt %>% inner_join(hits %>% transmute(Name = gene))

# Perform paired t-test on normalized expression SD between the two populations
t.test(genes_filt$AA_norm, genes_filt$EA_norm, paired=TRUE)
t.test(genes_hits$AA_norm, genes_hits$EA_norm, paired=TRUE)

# Define function to find expression SD in pop, a subset of AA individuals, for gene
pop_subset_sd = function(pop, gene) {
    gene_file = paste("data/fastqtl_sample_input/estimation/", pop, "/", gene, ".txt", sep="")
    if (!file.exists(gene_file)) {
        return(NA)
    }
    ind = read_tsv(gene_file, col_names=c("NWDID"))
    ind = ind %>% inner_join(samples)
    sd = reads %>% filter(Name == gene) %>% select(ind$TOR_ID) %>% sd
    return(sd)
}
pop_subset_sd_vect = Vectorize(pop_subset_sd)
genes_filt = genes_filt %>% mutate(Afr_sd = pop_subset_sd_vect("Afr", Name), het_sd = pop_subset_sd_vect("het", Name))
genes_filt %>% write.table("expression_std_dev.txt")

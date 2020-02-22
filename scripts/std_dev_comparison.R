library(tidyverse)
library(matrixStats)

metadata = read_csv("data/MESA_sample_info.csv")
samples = read_tsv("data/sample_participant_lookup.txt")
hits = read_tsv("results/v2/hits/hits_ascertainment_Eur.txt")
reads = read_tsv("data/mesa_TMM.expression.bed.gz")

# Filter for AA and EA samples by race information in individual metadata
AA = samples %>% left_join(metadata %>% select(NWDID, race1c) %>% unique) %>% filter(race1c == 3)
EA = samples %>% left_join(metadata %>% select(NWDID, race1c) %>% unique) %>% filter(race1c == 1)
# EA = EA %>% sample_n(AA %>% nrow) # Optionally downsample EA dataset until it has the same sample size as AA
AAreads = reads %>% select(AA$NWDID)
EAreads = reads %>% select(EA$NWDID)

# Within each population, compute standard deviation in gene reads for all genes
AAsd = AAreads %>% transmute(AA_sd = rowSds(as.matrix(.[AAreads %>% colnames])))
EAsd = EAreads %>% transmute(EA_sd = rowSds(as.matrix(.[EAreads %>% colnames])))

# Normalize expression SD measurement relative to AA expression SD
genes = bind_cols(reads %>% select('#chr', start, end, gene_id), AAsd, EAsd)
genes = genes %>% mutate(AA_norm = 1, EA_norm = EA_sd/AA_sd)

# Define function to find expression SD in pop, a subset of AA individuals, for gene
pop_subset_sd = function(pop, gene) {
    gene_file = paste("data/fastqtl_sample_input/estimation/", pop, "/", gene, ".txt", sep="")
    if (!file.exists(gene_file)) {
        return(NA)
    }
    ind = read_tsv(gene_file, col_names=c("NWDID"))
    ind = ind %>% inner_join(samples)
    sd = reads %>% filter(gene_id == gene) %>% select(ind$NWDID) %>% sd
    return(sd)
}
pop_subset_sd_vect = Vectorize(pop_subset_sd)
genes = genes %>% mutate(Afr_sd = pop_subset_sd_vect("Afr", gene_id), het_sd = pop_subset_sd_vect("het", gene_id), Eur_sd = pop_subset_sd_vect("Eur", gene_id))
genes %>% write.table("data/expression_std_dev.txt", row.names=FALSE, sep='\t')

# Run from home directory with at least 2G of memory

library(tidyverse)
library(qvalue)

dir.create("data/hits")
asc = read_tsv("data/fastqtl_output/merged_ascertainment_Eur.txt", col_names = c("gene", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "ref_factor", "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta")) # Heading pulled from Francois Aguet's FastQTL repo

# Apply Storey-Tibshirani FDR correction and MAF filter to ascertainment dataset to identify significant hits
asc = asc %>% mutate(fdr = qvalue(pval_beta)$qvalues)
hits = asc %>% filter(maf >=.05, fdr <= .10)
hits %>% write.table("data/hits/hits_ascertainment_Eur.txt", sep='\t', row.names=FALSE)

# Extract effect sizes (+ other info) for significant hits from estimation datasets
# Note that estimation datasets must be read and filtered sequentially otherwise you'll run out of memory
eur = read_tsv("data/fastqtl_output/merged_estimation_Eur.txt", col_names = c("gene", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"))
eur = eur %>% inner_join(hits %>% select(gene, variant))
eur %>% write.table("results/v2/hits/hits_estimation_Eur.txt", sep='\t', row.names=FALSE)

afr = read_tsv("data/fastqtl_output/merged_estimation_Afr.txt", col_names = c("gene", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"))
afr = afr %>% inner_join(hits %>% select(gene, variant))
afr %>% write.table("results/v2/hits/hits_estimation_Afr.txt", sep='\t', row.names=FALSE)

het = read_tsv("data/fastqtl_output/merged_estimation_het.txt", col_names = c("gene", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"))
het = het %>% inner_join(hits %>% select(gene, variant))
het %>% write.table("results/v2/hits/hits_estimation_het.txt", sep='\t', row.names=FALSE)

# Run with at least 30G of memory in home directory

library(tidyverse)

dir.create("data/hits")
asc = read_tsv("data/fastqtl_output/merged_ascertainment_Eur.txt", col_names = c("gene", "num_var", "beta_shape1", "beta_shape2", "true_df", "pval_true_df", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "ref_factor", "pval_nominal", "slope", "slope_se", "pval_perm", "pval_beta")) # Heading pulled from Francois Aguet's FastQTL repo
eur = read_tsv("data/fastqtl_output/merged_estimation_Eur.txt", col_names = c("gene", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"))
afr = read_tsv("data/fastqtl_output/merged_estimation_Afr.txt", col_names = c("gene", "variant", "tss_distance", "ma_samples", "ma_count", "maf", "pval_nominal", "slope", "slope_se"))

# Apply Bonferroni correction and MAF filter to ascertainment dataset to identify significant hits
total_tests = asc %>% pull("num_var") %>% sum
bonf_threshold = .05 / total_tests
hits = asc %>% filter(pval_beta < bonf_threshold, maf > .05)
hits %>% write.table("data/hits/hits_ascertainment_Eur.txt", sep='\t', row.names=FALSE)

# Extract effect sizes (+ other info) for significant hits from estimation datasets
eur = eur %>% inner_join(hits %>% select(gene, variant))
eur %>% write.table("data/hits/hits_estimation_Eur.txt", sep='\t', row.names=FALSE)
afr = afr %>% inner_join(hits %>% select(gene, variant))
afr %>% write.table("data/hits/hits_estimation_Afr.txt", sep='\t', row.names=FALSE)

library(tidyverse)
library(ggpubr)

setwd("~/sherlock/oak/eQTLs/results/likelihood_model/v4.2/")

# Specify column names of ANOVA df
anova_names = c("one_int", "race_int", "ganc_int", 
                "lanc_int", "one_beta", "anc_beta", 
                "delta", "anc_beta_AAEur", "delta_AAEur",
                "anc_beta_AAAfr", "delta_AAAfr",
                "anc_beta_EA", "delta_EA")

# Specify ordered file paths of different models
anova_paths = c("seq_center.exam/test.mode_fit_no_beta.ind_all",
                "race_Afr.race_Eur.seq_center.exam/test.mode_fit_no_beta.ind_all",
                "global_ancestry.genotype_PC1.race_Afr.race_Eur.seq_center.exam/test.mode_fit_no_beta.ind_all",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_fit_no_beta.ind_all",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_fit_single_beta.ind_all",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_no_interaction.ind_all",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_fit_delta.ind_all",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_no_interaction.ind_AAEur",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_fit_delta.ind_AAEur",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_no_interaction.ind_AAAfr",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_fit_delta.ind_AAAfr",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_no_interaction.ind_EA",
                "global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/test.mode_fit_delta.ind_EA")

anova_df = NULL
for (f in anova_paths) {
  tmp = read_tsv(paste0(f, ".optimize.anova_analysis.txt"))
  if (is.null(anova_df)) {
    anova_df = tmp
  } else {
    anova_df = inner_join(anova_df, tmp, by="gene")
  }
}
names(anova_df) = c("gene", anova_names)

anova_df %>% summarize(across(-gene, mean)) %>% data.frame

for (i in seq(2, 7)) {
  curr_mod = anova_names[i]
  prev_mod = anova_names[i - 1]
  diff_col = paste0("diff_", curr_mod)
  diff = (anova_df %>% pull(!! curr_mod)) - (anova_df %>% pull(!! prev_mod))
  anova_df = anova_df %>% mutate(!! diff_col := diff)
}
anova_df = anova_df %>% mutate(diff_delta_AAEur = delta_AAEur - anc_beta_AAEur)
anova_df = anova_df %>% mutate(diff_delta_AAAfr = delta_AAAfr - anc_beta_AAAfr)
anova_df = anova_df %>% mutate(diff_delta_EA = delta_EA - anc_beta_EA)
anova_df %>% select(starts_with("diff")) %>% summarize(across(everything(), mean)) %>% data.frame

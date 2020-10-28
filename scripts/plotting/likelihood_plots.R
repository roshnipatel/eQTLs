library(tidyverse)
library(ggpubr)
library(reshape2)

setwd("~/sherlock/oak/eQTLs/results/likelihood_model/v4.1/global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/")
dir.create("plots/", showWarnings = FALSE)

# Plot delta vs EM iteration
plot_delta_vs_iteration = function(delta_files) {
  delta_plots = list()
  i = 1
  for (f in delta_files) {
    df = read_tsv(f, col_names = c("Delta"))
    df = df %>% mutate(Iteration = seq(nrow(df)))
    p = ggplot(df, aes(Iteration, Delta)) + geom_point() + 
      ylim(0, 1) + theme_pubr() + scale_x_continuous(breaks=seq(0, 50, 10))
    delta_plots[[i]] = p
    i = i + 1
  }
  return(delta_plots)
}
pdf("plots/delta_vs_iteration_neg_control.pdf")
ggarrange(plotlist = plot_delta_vs_iteration(c("optimize_neg_control_delta.txt")), ncol=1, nrow=1)
dev.off()
pdf("plots/delta_vs_iteration_drop_asc.pdf")
ggarrange(plotlist = plot_delta_vs_iteration(c("optimize_drop_asc_delta.txt")), ncol=1, nrow=1)
dev.off()

# Plot likelihood vs delta
plot_likelihood = function(likelihood_file, bootstrap_file) {
  df = read_tsv(likelihood_file)
  min_val =  df %>% slice(which.max(Likelihood)) %>% pull(Delta)
  y_min = (df %>% summarise(min(Likelihood)))[[1]]
  y_max = (df %>% summarise(max(Likelihood)))[[1]]
  bootstrap_tbl = read.table(bootstrap_file)
  min_bootstrap = bootstrap_tbl[1,3]
  max_bootstrap = bootstrap_tbl[2,3]
  p = ggplot(df, aes(Delta, Likelihood)) + 
    annotate("rect", xmin=min_bootstrap, xmax=max_bootstrap, ymin=y_min - 20, ymax=y_max + 20, fill="#F5A76C") +
    geom_point() + 
    theme_pubr() + geom_vline(xintercept = min_val) + ylab("ln(Likelihood)") + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
    theme(text = element_text(size=20))
  return(p)
}
pdf("plots/likelihood_vs_delta_drop_asc.pdf")
plot_likelihood("likelihood_vs_delta_drop_asc.txt", "results_bootstrap_drop_asc_delta.summary.txt")
dev.off()
pdf("plots/likelihood_vs_delta_neg_control.pdf")
plot_likelihood("likelihood_vs_delta_neg_control.txt", "results_bootstrap_neg_control_delta.summary.txt")
dev.off()

# Plot distribution of bootstrap deltas
deltas = read_tsv("results_bootstrap_drop_asc_delta.all_values.txt", col_names = c("Delta"))
ggplot(deltas, aes(Delta)) + geom_histogram(binwidth = 0.005) + theme_pubr()

# Plot Afr vs Eur effects for EM run
source("~/sherlock/oak/eQTLs/scripts/plots/shared_fn.R")
new_betas = read_tsv("optimize_drop_asc_betas.txt")
new_betas = new_betas %>% transmute(Gene = gene, effect_Afr = beta_Afr,
                                    effect_Eur = beta_Eur)
pdf("plots/Afr_vs_Eur_effects.1dpca.pdf")
PCA_res = one_dim_PCA(new_betas %>% select(effect_Eur, effect_Afr))
ggplot(new_betas, aes(effect_Eur, effect_Afr)) + geom_point() + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], color='blue', size=1) + 
  geom_abline(slope=1, intercept=0) + 
  labs(x="European effect size", y="African effect size") + theme_pubr() + 
  xlim(-4, 7) + ylim(-4, 7) 
dev.off()
pdf("plots/Afr_vs_Eur_effects.linreg.pdf")
ggplot(new_betas, aes(effect_Eur, effect_Afr)) + geom_point() + 
  geom_abline(slope=1, intercept=0) + 
  geom_smooth(method='lm') +
  labs(x="European effect size", y="African effect size") + theme_pubr() + 
  xlim(-4, 7) + ylim(-4, 7) 
dev.off()
summary(lm(new_betas$effect_Afr ~ new_betas$effect_Eur))
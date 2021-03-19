library(tidyverse)
library(ggpubr)
library(reshape2)

setwd("~/sherlock/oak/eQTLs/results/likelihood_model/v4.2/global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/")
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
pdf("plots/likelihood_vs_delta.test.pdf")
plot_likelihood("test.likelihood_vs_delta.txt", "test.bootstrap_summary.txt")
dev.off()
pdf("plots/likelihood_vs_delta.control.pdf")
plot_likelihood("control.likelihood_vs_delta.txt", "control.bootstrap_summary.txt")
dev.off()

# Plot distribution of bootstrap deltas
deltas = read_tsv("results_bootstrap_drop_asc_delta.all_values.txt", col_names = c("Delta"))
ggplot(deltas, aes(Delta)) + geom_histogram(binwidth = 0.005) + theme_pubr()

# Plot Afr vs Eur effects for EM run
source("~/sherlock/oak/eQTLs/scripts/plotting/shared_fn.R")
new_betas = read_tsv("optimize_drop_asc_betas.txt")
new_betas = new_betas %>% transmute(gene, effect_Afr = beta_Afr,
                                    effect_Eur = beta_Eur)
PCA_res = one_dim_PCA(new_betas %>% select(effect_Eur, effect_Afr))
bootstrap_CI = quantile(bootstrap.TLS(new_betas %>% select(effect_Eur, effect_Afr)), c(0.025, 0.975))
plot_df = new_betas %>% mutate(dummy_x = seq(from=-6, to=9, length.out=nrow(new_betas)), 
                            CI_min = dummy_x * bootstrap_CI[[1]] + PCA_res[[2]], 
                            CI_max = dummy_x * bootstrap_CI[[2]] + PCA_res[[2]])
pdf("plots/Afr_vs_Eur_effect_sizes.TLS.pdf")
ggplot(plot_df) + 
  geom_point(mapping=aes(effect_Eur, effect_Afr)) + 
  labs(x="European effect size", y="African effect size") + theme_pubr() +
  geom_ribbon(aes(x=dummy_x, ymin=CI_max, ymax=CI_min), fill="#F5A76C") +
  geom_abline(slope=1, intercept=0, linetype="dotted") + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], linetype="solid") + 
  coord_cartesian(xlim=c(-5, 7), ylim=c(-5, 7)) + 
  theme(text = element_text(size=20))
dev.off()

# Compare betas from full model vs initial estimation of anc-specific effect sizes
merged = set_up("~/sherlock/oak/eQTLs/results/QTL_calling/v4.0/hits/") %>% 
  select(gene, effect_Afr, effect_Eur) %>% pivot_longer(-gene)
new_betas = new_betas %>% pivot_longer(-gene)
merged = merged %>% inner_join(new_betas, by=c("gene", "name"), suffix = c("_OLS", "_model"))
p_Eur = ggplot(merged %>% filter(name == "effect_Eur"), aes(value_OLS, value_model)) +
  geom_point() + theme_pubr() + theme(text = element_text(size=20)) +
  geom_abline(slope = 1, intercept = 0, linetype="dotted") + xlab("OLS estimate") + 
  ylab("likelihood model estimate") + ggtitle("European effects")
p_Afr = ggplot(merged %>% filter(name == "effect_Afr"), aes(value_OLS, value_model)) +
  geom_point() + theme_pubr() + theme(text = element_text(size=20)) +
  geom_abline(slope = 1, intercept = 0, linetype="dotted") + xlab("OLS estimate") + 
  ylab("likelihood model estimate") + ggtitle("African effects")
pdf("plots/compare_model_betas.pdf")
ggarrange(p_Eur, p_Afr)
dev.off()

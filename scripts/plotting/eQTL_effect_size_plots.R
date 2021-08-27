library(tidyverse)
library(ggpubr)
source("~/sherlock/oak/admixedPGS/eQTLs/scripts/plotting/shared_fn.R")

# Set working directory to local directory containing paper figures
setwd("Google Drive/stanford/admixed PGS/figures/")
dir.create("plots", showWarnings = FALSE)

############################## For MESA ##############################
merged = inner_join(read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/QTL_calling/v4.4/hits/hits_reestimation_primary_AA.txt"), 
                    read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/QTL_calling/v4.4/hits/hits_reestimation_primary_EA.txt"),
                    by=c("gene", "ID"), suffix=c("_AA", "_EA"))
merged = inner_join(read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/QTL_calling/v4.4/hits/hits_ascertainment_Eur.txt"),
                    merged, by=c("gene", "ID"))

# Plot AA vs EA effect sizes and total least squares fit
PCA_res = one_dim_PCA(merged %>% select(effect_EA, effect_AA))
bootstrap_CI = quantile(bootstrap.TLS(merged %>% select(effect_EA, effect_AA)), c(0.025, 0.975))
plot_df = merged %>% mutate(dummy_x = seq(from=-6, to=9, length.out=nrow(merged)), 
                            CI_min = dummy_x * bootstrap_CI[[1]] + PCA_res[[2]], 
                            CI_max = dummy_x * bootstrap_CI[[2]] + PCA_res[[2]])
slope_label = paste0("Slope = ", round(PCA_res[[1]], digits=2))
r2_label = paste0("R2 = ", round(cor(merged$effect_EA, merged$effect_AA), digits=2))
ggplot(plot_df) + 
  labs(x="European-American effect size", y="African-American effect size") + theme_pubr() +
  geom_ribbon(aes(x=dummy_x, ymin=CI_max, ymax=CI_min), fill="#F5A76C") +
  geom_abline(slope=1, intercept=0, linetype="dashed") + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], linetype="solid") + 
  geom_point(mapping=aes(effect_EA, effect_AA)) + 
  coord_cartesian(xlim=c(-5, 7.5), ylim=c(-5, 7.5)) + 
  ggtitle("a) Gene expression") +
  theme(axis.text = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5)) + 
  annotate(geom="text", x=-3.5, y=7.5, label=plot_label, size=10) 
ggsave("plots/AA_vs_EA_effect_sizes.TLS.pdf", height=10, width=10, units="in")

# Plot AA vs EA effect sizes and regress with OLS and inverse variance weighting
weights = 1/sqrt((merged %>% pull(se_AA))^2 + (merged %>% pull(se_EA))^2) # Define weights for inverse variance weighting
mod = lm(effect_AA ~ effect_EA, merged)
ggplot(merged, aes(effect_EA, effect_AA)) + geom_point() + 
  geom_abline(intercept=0, slope=1, linetype="dashed") + 
  stat_smooth(method=lm, color='black', aes(weight=weights), 
              fill="#F5A76C", alpha = 1, fullrange = TRUE, size=0.5) + 
  labs(x="European-American effect size", y="African-American effect size") + 
  xlim(-5, 7) + ylim(-5, 7) + theme_pubr() +
  ggtitle("a) Gene expression") +
  theme(axis.text = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5))
ggsave("plots/AA_vs_EA_effect_sizes.OLS.pdf", height=10, width=10, units="in")

############################## For MVP ##############################
rm(merged)
merged = read_delim("mvp_data/betaA_betaE_maf0.001_v6_122_snps.txt", delim=" ")

# Plot AA vs EA effect sizes and total least squares fit
PCA_res = one_dim_PCA(merged %>% select(effect_EA, effect_AA))
bootstrap_CI = quantile(bootstrap.TLS(merged %>% select(effect_EA, effect_AA)), c(0.025, 0.975))
plot_df = merged %>% mutate(dummy_x = seq(from=-6, to=9, length.out=nrow(merged)), 
                            CI_min = dummy_x * bootstrap_CI[[1]] + PCA_res[[2]], 
                            CI_max = dummy_x * bootstrap_CI[[2]] + PCA_res[[2]])
ggplot(plot_df) + 
  labs(x="European-American effect size", y="African-American effect size") + theme_pubr() +
  geom_ribbon(aes(x=dummy_x, ymin=CI_max, ymax=CI_min), fill="#F5A76C") +
  geom_abline(slope=1, intercept=0, linetype="dashed") + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], linetype="solid") + 
  geom_point(mapping=aes(effect_EA, effect_AA)) + 
  coord_cartesian(xlim=c(-0.15, 0.21), ylim=c(-0.15, 0.21)) + 
  ggtitle("b) LDL-C") +
  theme(axis.text = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5))
ggsave("plots/LDLC.AA_vs_EA_effect_sizes.TLS.pdf", height=10, width=10, units="in")

# Plot AA vs EA effect sizes and regress with OLS and inverse variance weighting
weights = 1/sqrt((merged %>% pull(se_AA))^2 + (merged %>% pull(se_EA))^2) # Define weights for inverse variance weighting
mod = lm(effect_AA ~ effect_EA, merged)
ggplot(merged, aes(effect_EA, effect_AA)) + 
  geom_abline(intercept=0, slope=1, linetype="dashed") + 
  stat_smooth(method=lm, color='black', aes(weight=weights), 
              fill="#F5A76C", alpha = 1, fullrange = TRUE, size=0.5) + 
  labs(x="European-American effect size", y="African-American effect size") + 
  coord_cartesian(xlim=c(-0.15, 0.21), ylim=c(-0.15, 0.21)) + theme_pubr() +
  ggtitle("b) LDL-C") + geom_point() + 
  theme(axis.text = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5))
ggsave("plots/LDLC.AA_vs_EA_effect_sizes.OLS.pdf", height=10, width=10, units="in")

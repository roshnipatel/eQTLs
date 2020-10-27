library(tidyverse)
library(ggpubr)

# Specify working directory
setwd("~/sherlock/oak/eQTLs/results/QTL_calling/v4.0/")
dir.create("plots", showWarnings = FALSE)

source("~/sherlock/oak/eQTLs/scripts/plotting/shared_fn.R")

merged = set_up("hits/")
# merged = merged %>% mutate(effect_Eur = effect_val, se_Eur = se_val) # Uncomment if using Afr ascertainment (v4.1)
merged = merged %>% mutate(t_Eur = effect_Eur / se_Eur, t_Afr = effect_Afr / se_Afr)


# Plot QQ plots
QQ_plot = function(eQTL_path, plot_title) {
  tib = read_tsv(eQTL_path)
  tib = tib %>% arrange(pval) %>% mutate(expected = ppoints(pval))
  p = ggplot(tib, aes(-log10(expected), -log10(pval))) + geom_point() + 
    geom_abline(aes(intercept=0, slope=1)) + theme_pubr() +
    labs(title = plot_title)
  return(p)
}

# Obtain bootstrap confidence intervals for total least squares
bootstrap.TLS <- function(tib) {
  n = nrow(tib)
  slopes = rep(0, n)
  for (i in 1:n) {
    samp = tib %>% sample_n(n, TRUE)
    samp_slope = one_dim_PCA(samp)[[1]]
    slopes[i] = samp_slope
  }
  return(slopes)
}

# Plot Afr vs Eur t-statistics and total least squares fit
PCA_res = one_dim_PCA(merged %>% select(t_Eur, t_Afr))
bootstrap_CI = quantile(bootstrap.TLS(merged %>% select(t_Eur, t_Afr)), c(0.025, 0.975))
plot_df = merged %>% mutate(dummy_x = seq(from=-30, to=30, length.out=nrow(merged)), 
                            CI_min = dummy_x * bootstrap_CI[[1]] + PCA_res[[2]], 
                            CI_max = dummy_x * bootstrap_CI[[2]] + PCA_res[[2]])
line_df = data.frame(s = c(PCA_res[[1]], 1), ic = c(PCA_res[[2]], 0), legend = c("total least squares fit", "y = x"))
line_df$legend = factor(line_df$legend, levels=c("total least squares fit", "y = x"))
ggplot(plot_df) + 
  labs(x="European t-statistic", y="African t-statistic") + theme_pubr() +
  theme(legend.position = "right") +
  geom_ribbon(aes(x=dummy_x, ymin=CI_max, ymax=CI_min), fill="#F5A76C") +
  geom_abline(data=line_df, mapping=aes(slope=s, intercept=ic, linetype=legend)) + 
  geom_point(aes(t_Eur, t_Afr)) + 
  coord_cartesian(xlim=c(-25, 25), ylim=c(-25, 25)) + 
  theme(text = element_text(size=20))
ggsave("plots/Afr_vs_Eur_tstats.both.pdf")

# Plot Afr vs Eur effect sizes and total least squares fit
PCA_res = one_dim_PCA(merged %>% select(effect_Eur, effect_Afr))
bootstrap_CI = quantile(bootstrap.TLS(merged %>% select(effect_Eur, effect_Afr)), c(0.025, 0.975))
plot_df = merged %>% mutate(dummy_x = seq(from=-6, to=9, length.out=nrow(merged)), 
                            CI_min = dummy_x * bootstrap_CI[[1]] + PCA_res[[2]], 
                            CI_max = dummy_x * bootstrap_CI[[2]] + PCA_res[[2]])
line_df = data.frame(s = c(PCA_res[[1]], 1), ic = c(PCA_res[[2]], 0), legend = c("total least squares fit", "y = x"))
line_df$legend = factor(line_df$legend, levels=c("total least squares fit", "y = x"))
ggplot(plot_df) + 
  labs(x="European effect size", y="African effect size") + theme_pubr() +
  theme(legend.position = "right") +
  geom_ribbon(aes(x=dummy_x, ymin=CI_max, ymax=CI_min), fill="#F5A76C") +
  geom_abline(data=line_df, mapping=aes(slope=s, intercept=ic, linetype=legend)) + 
  geom_point(aes(effect_Eur, effect_Afr)) + 
  coord_cartesian(xlim=c(-5, 7), ylim=c(-5, 7)) + 
  theme(text = element_text(size=20))
ggsave("plots/Afr_vs_Eur_effect_sizes.both.pdf")

ggplot(plot_df) + geom_ribbon(aes(x=dummy_x, ymin=CI_max, ymax=CI_min), fill="#F5A76C")

# Plot Afr vs Eur effect sizes and regress with OLS and inverse variance weighting
weights = 1/sqrt((merged %>% pull(se_Afr))^2 + (merged %>% pull(se_Eur))^2) # Define weights for inverse variance weighting
pdf("plots/Afr_vs_Eur_effect_sizes.OLS.pdf")
p_AfrvsEur = ggplot(merged, aes(effect_Eur, effect_Afr)) + geom_point() + 
  geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1)) + 
  labs(x="European effect size", y="African effect size") + theme_pubr() + 
  stat_cor(method = "pearson", label.x = -2, label.y = 3.5, size=6) +
  stat_regline_equation(label.x=-2, label.y = 4, size=6) + 
  xlim(-2, 5) + ylim(-2, 5)
p_EurvsAfr = ggplot(merged, aes(effect_Afr, effect_Eur)) + geom_point() + 
  geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1)) + 
  labs(x="African effect size", y="European effect size") + theme_pubr() + 
  stat_cor(method = "pearson", label.x = -2, label.y = 3.5, size=6) +
  stat_regline_equation(label.x=-2, label.y = 4, size=6) + 
  xlim(-2, 5) + ylim(-2, 5)
ggarrange(p_AfrvsEur, p_EurvsAfr)
dev.off()

# Plot Afr vs Eur effect sizes and regress with 1d PCA
pdf("plots/Afr_vs_Eur_effect_sizes.1dPCA.pdf")
PCA_res = one_dim_PCA(merged %>% select(effect_Eur, effect_Afr))
ggplot(merged, aes(effect_Eur, effect_Afr)) + geom_point() + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], color='blue', size=1) + 
  geom_abline(slope=1, intercept=0) + 
  labs(x="European effect size", y="African effect size") + theme_pubr() +
  xlim(-2, 5) + ylim(-2, 5)
dev.off()

# Plot hypothesis testing and regress with OLS
pdf("plots/effect_dif_wtavg.OLS.pdf")
ggplot(merged, aes(exp_dif_wtavg, obs_dif_het)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + 
  labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + 
  geom_smooth(method=lm, color='blue') + theme_pubr() + 
  stat_cor(method = "pearson", label.x = -3, label.y = .75, size=6) +
  stat_regline_equation(label.x=-3, label.y = 1, size=6) +
  xlim(-3, 1.5) + ylim(-3, 1.5)
ggplot(merged, aes(obs_dif_het, exp_dif_wtavg)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + 
  labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + 
  geom_smooth(method=lm, color='blue') + theme_pubr() + 
  stat_cor(method = "pearson", label.x = -3, label.y = .75, size=6) +
  stat_regline_equation(label.x=-3, label.y = 1, size=6) +
  xlim(-3, 1.5) + ylim(-3, 1.5)
dev.off()

# Plot hypothesis testing and regress with 1d PCA
pdf("plots/effect_dif_wtavg.1dPCA.pdf")
PCA_res = one_dim_PCA(merged %>% select(exp_dif_wtavg, obs_dif_het))
ggplot(merged, aes(exp_dif_wtavg, obs_dif_het)) + geom_point() + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], color='blue', size=1) + 
  geom_abline(slope=1, intercept=0) + 
  labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + 
  theme_pubr() + xlim(-3, 1.5) + ylim(-3, 1.5)
dev.off()

# Plot negative control and regress with OLS
pdf("plots/validation_vs_wtavg.OLS.pdf")
ggplot(merged, aes(exp_dif_wtavg, obs_dif_val)) + geom_point() + 
  geom_abline(aes(slope=1, intercept=0)) + 
  labs(x="Exp. het. dif. from Eur. effect size", y="Dif. between Eur. effect size in validation and estimation sets") + 
  geom_smooth(method=lm, color='blue') + theme_pubr() +
  stat_cor(method = "pearson", label.x = -3, label.y = 1.05, size=6) +
  stat_regline_equation(label.x=-3, label.y = 1.3, size=6) + 
  xlim(-3, 1.5) + ylim(-3, 1.5)
dev.off()

# Plot negative control and regress with 1d PCA
PCA_res = one_dim_PCA(merged %>% select(exp_dif_wtavg, obs_dif_val))
pdf("plots/validation_vs_wtavg.OLS.pdf")
ggplot(merged, aes(exp_dif_wtavg, obs_dif_val)) + geom_point() + 
  geom_abline(slope=PCA_res[[1]], intercept=PCA_res[[2]], color='blue', size=1) + 
  geom_abline(slope=1, intercept=0) + 
  labs(x="Exp. dif. between het effect size and Eur. effect size", y="Obs. dif. between validation and estimation Eur. effect size") + 
  theme_pubr() + xlim(-3, 1.5) + ylim(-3, 1.5)
dev.off()

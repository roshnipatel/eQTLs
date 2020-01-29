# Run from home directory
# IMPORTANT: change file paths

library(tidyverse)

het = read_tsv("results/v2/hits/hits_estimation_het.txt")
afr = read_tsv("results/v2/hits/hits_estimation_Afr.txt")
eur = read_tsv("results/v2/hits/hits_estimation_Eur.txt")

merged = het %>% inner_join(afr, by=c("gene", "variant"), suffix=c("_het", "_Afr")) %>% inner_join(eur, by=c("gene", "variant"), suffix=c("", "_Eur")) # Note that the last suffix command doesn't work properly so Eur columns have no suffix appended
merged = merged %>% mutate(AfrEur_avg = (slope_Afr + slope) / 2)
merged = merged %>% mutate(varAfr = maf_Afr * (1 - maf_Afr), varEur = maf * (1 - maf))
merged = merged %>% mutate(alpha = varAfr / (varAfr + varEur))
merged = merged %>% mutate(AfrEur_wtavg = alpha * slope_Afr + (1 - alpha) * slope)
merged = merged %>% mutate(obs_dif = slope_het - slope, exp_dif_avg = AfrEur_avg - slope, exp_dif_wtavg = AfrEur_wtavg - slope)

# Plot Afr effect sizes vs Eur effect sizes without confidence intervals
tib = merged %>% filter(maf_Afr > .05, maf > .05) # Filter for eQTLs with MAF of at least .05 in both estimation datasets
weights = 1/sqrt((tib %>% pull(slope_se_Afr))^2 + (tib %>% pull(slope_se))^2) # Define weights for inverse variance weighting
pdf("results/v2/plots/Afr_vs_Eur_effect_sizes.pdf")
ggplot(tib, aes(slope, slope_Afr)) + geom_point() + geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1)) + labs(x="European effect size", y="African effect size")
dev.off()

# Plot Afr effect sizes vs Eur effect sizes with confidence intervals
tib = merged %>% filter(maf_Afr > .05, maf > .05) # Filter for eQTLs with MAF of at least .05 in both estimation datasets
weights = 1/sqrt((tib %>% pull(slope_se_Afr))^2 + (tib %>% pull(slope_se))^2) # Define weights for inverse variance weighting
pdf("results/v2/plots/Afr_vs_Eur_effect_sizes_CI.pdf")
ggplot(tib, aes(slope, slope_Afr)) + geom_segment(tib, mapping = aes(x=slope -1.96*slope_se, xend=slope + 1.96 * slope_se, y=slope_Afr, yend = slope_Afr), color='gray', alpha=.9, size=.8) + geom_segment(tib, mapping = aes(x=slope, xend=slope, y=slope_Afr - 1.96 * slope_se_Afr, yend = slope_Afr + 1.96 * slope_se_Afr), color='gray', alpha=.9, size=.8) + geom_point() + geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1)) + labs(x="European effect size", y="African effect size")
dev.off()

# Plot observed vs expected deviation from European effect sizes in ancestry-heterozygous individuals
# Inverse variance weighting assumes no covariance between effect sizes in difference ancestry backgrounds

# Use average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
weights = merged %>% transmute(varY = slope_se_het^2 + slope_se^2, varX = .25 * (slope_se_Afr ^2 + slope_se ^2)) %>% mutate(weight = 1/sqrt(varX + varY))
pdf("results/v2/plots/slope_dif_avg_invvar_weight.pdf")
ggplot(merged, aes(exp_dif_avg, obs_dif)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + geom_smooth(method=lm, color='blue', aes(weight = weights$weight))
dev.off()

# Use allele frequency-weighted average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
weights = merged %>% transmute(varY = slope_se_het^2 + slope_se^2, varX = (varAfr^2 * (slope_se_Afr ^2) + varEur^2 * (slope_se ^2)) / (varAfr + varEur)^2 ) %>% mutate(weight = 1/sqrt(varX + varY))
pdf("results/v2/plots/slope_dif_wtavg_invvar_weight.pdf")
ggplot(merged, aes(exp_dif_wtavg, obs_dif)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + geom_smooth(method=lm, color='blue', aes(weight=weights$weight))
dev.off()

# Plot observed vs expected deviation from European effect sizes in ancestry-heterozygous individuals
# Do NOT use weights on regression

# Use average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
pdf("results/v2/plots/slope_dif_avg_no_weight.pdf")
ggplot(merged, aes(exp_dif_avg, obs_dif)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + geom_smooth(method=lm, color='blue')
dev.off()

# Use allele frequency-weighted average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
pdf("results/v2/plots/slope_dif_wtavg_no_weight.pdf")
ggplot(merged, aes(exp_dif_wtavg, obs_dif)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + geom_smooth(method=lm, color='blue')
dev.off()

# Plot absolute value of observed vs expected deviation from European effect sizes in ancestry-heterozygous individuals
# Use average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
pdf("results/v2/plots/slope_dif_avg_abs.pdf")
ggplot(merged, aes(abs(exp_dif_avg), abs(obs_dif))) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") 
dev.off()

# Use allele frequency-weighted average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
pdf("results/v2/plots/slope_dif_wtavg_abs.pdf")
ggplot(merged, aes(abs(exp_dif_wtavg), abs(obs_dif))) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size")
dev.off()

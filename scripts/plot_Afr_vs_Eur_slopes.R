# Run from home directory

library(tidyverse)

anchet = read_tsv("data/hits/hits_estimation_anc_het.txt")
afr = read_tsv("data/hits/hits_estimation_Afr.txt")
eur = read_tsv("data/hits/hits_estimation_Eur.txt")

merged = anchet %>% inner_join(afr, by=c("gene", "variant"), suffix=c("_anc_het", "_Afr")) %>% inner_join(eur, by=c("gene", "variant"), suffix=c("", "_Eur")) # Note that the last suffix command doesn't work properly so Eur columns have no suffix appended
merged = merged %>% mutate(AfrEur_avg = (slope_Afr + slope) / 2)
merged = merged %>% mutate(AfrEur_wtavg = (maf_Afr / (maf_Afr + maf)) * slope_Afr + (maf / (maf_Afr + maf)) * slope)
merged = merged %>% mutate(obs_dif = slope_anc_het - slope, exp_dif_avg = AfrEur_avg - slope, exp_dif_wtavg = AfrEur_wtavg - slope)

# Plot Afr effect sizes vs Eur effect sizes without confidence intervals
tib = merged %>% filter(maf_Afr > .05, maf > .05) # Filter for eQTLs with MAF of at least .05 in both estimation datasets
weights = 1/sqrt((tib %>% pull(slope_se_Afr))^2 + (tib %>% pull(slope_se))^2) # Define weights for inverse variance weighting
pdf("plots/Afr_vs_Eur_effect_sizes.pdf")
ggplot(tib, aes(slope, slope_Afr)) + geom_point() + geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1)) + labs(x="European effect size", y="African effect size")
dev.off()

# Plot Afr effect sizes vs Eur effect sizes with confidence intervals
tib = merged %>% filter(maf_Afr > .05, maf > .05) # Filter for eQTLs with MAF of at least .05 in both estimation datasets
weights = 1/sqrt((tib %>% pull(slope_se_Afr))^2 + (tib %>% pull(slope_se))^2) # Define weights for inverse variance weighting
pdf("plots/Afr_vs_Eur_effect_sizes_CI.pdf")
ggplot(tib, aes(slope, slope_Afr)) + geom_segment(tib, mapping = aes(x=slope -1.96*slope_se, xend=slope + 1.96 * slope_se, y=slope_Afr, yend = slope_Afr), color='gray', alpha=.9, size=.8) + geom_segment(tib, mapping = aes(x=slope, xend=slope, y=slope_Afr - 1.96 * slope_se_Afr, yend = slope_Afr + 1.96 * slope_se_Afr), color='gray', alpha=.9, size=.8) + geom_point() + geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1)) + labs(x="European effect size", y="African effect size")
dev.off()

# Plot observed vs expected deviation from European effect sizes in ancestry-heterozygous individuals

# Use average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
weights = merged %>% transmute(varX = slope_se_anc_het^2 + slope_se^2, varY = .25 * (slope_se_Afr ^2 + slope_se ^2)) %>% mutate(weight = 1/sqrt(varX + varY))
pdf("plots/slope_dif_avg.pdf")
ggplot(merged, aes(exp_dif_avg, obs_dif)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + geom_smooth(method=lm, color='blue', weight = weights$weight)
dev.off()

# Use allele frequency-weighted average of Afr and Eur effect sizes to represent expected effect size in ancestry-heterozygous
weights = merged %>% transmute(varX = slope_se_anc_het^2 + slope_se^2, varY = (maf_Afr / (maf_Afr + maf))^2 * (slope_se_Afr ^2) + (maf / (maf_Afr + maf))^2 * (slope_se ^2)) %>% mutate(weight = 1/sqrt(varX + varY)
pdf("plots/slope_dif_wtavg.pdf")
ggplot(merged, aes(exp_dif_wtavg, obs_dif)) + geom_point() + geom_abline(aes(slope=1, intercept=0)) + labs(x="Exp. dif. from Eur. effect size", y="Obs. dif. from Eur. effect size") + geom_smooth(method=lm, color='blue', weight=weights$weight)
dev.off()

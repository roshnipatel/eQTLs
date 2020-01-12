# Run from home directory

library(tidyverse)

afr = read_tsv("data/hits/hits_estimation_Afr.txt")
eur = read_tsv("data/hits/hits_estimation_Eur.txt")

tib = inner_join(afr, eur, by=c("gene", "variant"), suffix=c("_Afr", "_Eur"))
tib = tib %>% filter(maf_Afr > .05, maf_Eur > .05) # Filter for eQTLs with MAF of at least .05 in both estimation datasets

weights = 1/sqrt((tib %>% pull(slope_se_Afr))^2 + (tib %>% pull(slope_se_Eur))^2) # Define weights for inverse variance weighting

pdf("plots/Afr_vs_Eur_effect_sizes.pdf")
ggplot(tib, aes(slope_Eur, slope_Afr)) + geom_segment(tib, mapping = aes(x=slope_Eur -1.96*slope_se_Eur, xend=slope_Eur + 1.96 * slope_se_Eur, y=slope_Afr, yend = slope_Afr), color='gray', alpha=.9, size=.8) + geom_segment(tib, mapping = aes(x=slope_Eur, xend=slope_Eur, y=slope_Afr - 1.96 * slope_se_Afr, yend = slope_Afr + 1.96 * slope_se_Afr), color='gray', alpha=.9, size=.8) + geom_point() + geom_smooth(method=lm, color='blue', aes(weight=weights)) + geom_abline(aes(intercept=0, slope=1))
dev.off()

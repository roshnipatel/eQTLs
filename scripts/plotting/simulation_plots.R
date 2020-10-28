library(tidyverse)
library(ggpubr)

setwd("~/sherlock/oak/eQTLs/results/")
dir.create("simulations/plots", showWarnings = FALSE)

sim_data = read_tsv("likelihood_model/v4.0/global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1/")
sim_betas = sim_data %>% select(gene, beta_Afr, beta_Eur, effect_Asc, effect_Afr, effect_Eur) %>% unique()

# Plot density of simulated true African and European betas
tib = sim_betas %>% select(beta_Afr, beta_Eur) %>% pivot_longer(c(beta_Afr, beta_Eur))
pdf("simulations/plots/true_beta_distribution.pdf")
ggplot(tib, aes(value)) + geom_density(aes(color=name)) + theme_pubr()
dev.off()

# Plot density of simulated and empirical estimated betas
real_asc = read_tsv("QTL_calling/v4.0/hits/hits_ascertainment_Eur.txt") %>% select(effect) %>% mutate(name="empirical")
tib = real_asc %>% bind_rows(sim_betas %>% select(effect_Asc) %>% transmute(name = "simulated", effect = effect_Asc))
pdf("simulations/plots/sim_vs_empirical_ascertainment.pdf")
ggplot(tib, aes(effect)) + geom_density(aes(color=name)) + theme_pubr()
dev.off()

real_eur = read_tsv("QTL_calling/v4.0/hits/hits_reestimation_primary_Eur.txt")
tib = real_eur %>% select(effect) %>% mutate(name="empirical") %>% bind_rows(sim_betas %>% select(effect_Eur) %>% transmute(name = "simulated", effect = effect_Eur))
pdf("simulations/plots/sim_vs_empirical_Eur_estimation.pdf")
ggplot(tib, aes(effect)) + geom_density(aes(color=name)) + theme_pubr()
dev.off()

real_afr = read_tsv("QTL_calling/v4.0/hits/hits_reestimation_primary_Afr.txt") 
tib = real_afr %>% select(effect) %>% mutate(name="empirical") %>% bind_rows(sim_betas %>% select(effect_Afr) %>% transmute(name = "simulated", effect = effect_Afr))
pdf("simulations/plots/sim_vs_empirical_Afr_estimation.pdf")
ggplot(tib, aes(effect)) + geom_density(aes(color=name)) + theme_pubr()
dev.off()

# Plot distribution of simulated global ancestry in an ADMIXTURE-type plot
n_idv = 150
sim_ganc = tibble(Ancestry = rbeta(n = n_idv, shape1 = 7.9, shape2 = 2.1)) %>% arrange(desc(Ancestry)) %>% mutate(Individual = seq(n_idv))
ggplot(sim_ganc, aes(Individual, Ancestry)) + geom_point() + theme_pubr() + ylim(0, 1)

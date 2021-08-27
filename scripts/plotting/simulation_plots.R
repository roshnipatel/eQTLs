library(tidyverse)
library(ggpubr)

# Set working directory to local directory containing paper figures
setwd("Google Drive/stanford/admixed PGS/figures/")
dir.create("plots", showWarnings = FALSE)

# Plot values of delta
sim_deltas = read_tsv("~/sherlock/scratch/admixedPGS/eQTLs/data/model_fitting/simulated_optimization/vA0.08.vE0.15.errA0.17.errE0.12/combined_optimization_outputs.txt")
ggplot(sim_deltas, aes(factor(true_delta), estimated_delta)) + 
  geom_violin() + theme_pubr() + 
  theme(axis.text.x = element_text(size=25), 
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5)) +
  labs(x="Simulated Delta", y="Estimated Delta") +
  ggtitle("a)")
ggsave("plots/estimated_vs_simulated_deltas.pdf", height=10, width=10, units="in")

# Read simulated data
sim_betas = read_tsv("~/sherlock/scratch/admixedPGS/eQTLs/data/model_fitting/simulated_data/vA0.08.vE0.15.errA0.17.errE0.12/betas.txt")

# Plot density of simulated and empirical estimated betas
real_asc = read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/QTL_calling/v4.4/hits/hits_ascertainment_Eur.txt")
tib = real_asc %>% select(effect) %>% mutate(name="empirical") %>% bind_rows(sim_betas %>% select(effect_Asc) %>% transmute(name = "simulated", effect = effect_Asc))
asc_plot = ggplot(tib, aes(effect)) + geom_density(aes(color=name)) + 
  theme_pubr() + ggtitle("European-American ascertainment") +
  theme(axis.text.x = element_text(size=25), 
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5)) +
  ggtitle("b)")
ggsave("plots/asc.estimated_vs_simulated_beta_distributions.pdf", height=10, width=10, units="in")

real_eur = read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/QTL_calling/v4.4/hits/hits_reestimation_primary_Eur.txt")
tib = real_eur %>% select(effect) %>% mutate(name="empirical") %>% bind_rows(sim_betas %>% select(effect_Eur) %>% transmute(name = "simulated", effect = effect_Eur))
eur_plot = ggplot(tib, aes(effect)) + geom_density(aes(color=name)) + 
  theme_pubr() + ggtitle("European-American re-estimation") + 
  theme(axis.text.x = element_text(size=25), 
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5)) +
  ggtitle("c)")
ggsave("plots/eur.estimated_vs_simulated_beta_distributions.pdf", height=10, width=10, units="in")

real_afr = read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/QTL_calling/v4.4/hits/hits_reestimation_primary_Afr.txt") 
tib = real_afr %>% select(effect) %>% mutate(name="empirical") %>% bind_rows(sim_betas %>% select(effect_Afr) %>% transmute(name = "simulated", effect = effect_Afr))
afr_plot = ggplot(tib, aes(effect)) + geom_density(aes(color=name)) + 
  theme_pubr() + ggtitle("African-American re-estimation") +
  theme(axis.text.x = element_text(size=25), 
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5)) +
  ggtitle("d)")
ggsave("plots/afr.estimated_vs_simulated_beta_distributions.pdf", height=10, width=10, units="in")

# Plot distribution of simulated global ancestry in an ADMIXTURE-type plot
n_idv = 319
# sim_ganc = tibble(Ancestry = rbeta(n = n_idv, shape1 = 7.9, shape2 = 2.1)) %>% arrange(desc(Ancestry)) %>% mutate(Individual = seq(n_idv))
sim_ganc = tibble(Ancestry = rbeta(n = n_idv, shape1 = 7.9, shape2 = 2.1))
true_ganc = read_tsv("~/sherlock/oak/admixedPGS/LocalAncestry/results/v2/combined_global_anc_frac.txt")
ggplot(sim_ganc, aes(Ancestry)) + geom_histogram() + theme_pubr()

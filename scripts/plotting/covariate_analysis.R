library(tidyverse)
library(ggpubr)

setwd("~/sherlock/oak/eQTLs/results/covariate_analysis/v2/")

# Define directory in which to store plots
dir.create("plots/", showWarnings = FALSE)

# Read sorted eigenvalue file and plot scree plot
plot_scree = function(filepath, n_eig) {
  eig = read_tsv(filepath, col_names =c("Eigenvalue"))
  eig = eig %>% mutate(Component = seq(1, nrow(eig)))
  p = ggplot(eig %>% head(n_eig), aes(Component, Eigenvalue)) + geom_point() + theme_pubr()
  return(p)
}
pdf("plots/genotype_PC.top30_scree_plot.pdf")
plot_scree("mesa.prune.maf05.vcf.eigenvalues.txt", 30)
dev.off()
pdf("plots/expression_PC.top30_scree_plot.pdf")
plot_scree("mesa.expression.bed.eigenvalues.txt", 30)
dev.off()

# Read covariate file
pc_tib = read_tsv("mesa.sample_covariates_deidentified.txt")

# Plot 2 PCs and color by a specified covariate column
plot_pc = function(tib, pc_x, pc_y, fill_variable) {
  p = ggplot(tib, aes_string(pc_x, pc_y)) + geom_point(aes_string(color = fill_variable)) + theme_pubr()
  return(p)
}
plot_pc(pc_tib, "genotype_PC0", "genotype_PC1", "race")
plot_pc(pc_tib, "expression_PC0", "expression_PC1", "exam")
plot_pc(pc_tib, "expression_PC1", "expression_PC2", "seq_center")

# Plot correlation matrix for covariates
plot_heatmap = function(tib) {
  cor_tib = as_tibble(tib %>% cor)
  long_cor_tib = cor_tib %>% mutate(covariate = colnames(cor_tib)) %>% pivot_longer(-covariate) %>% 
    mutate(value = abs(value))
  p = ggplot(long_cor_tib, aes(name, covariate, fill=value)) + scale_fill_distiller(palette = "RdPu") + 
    xlab("covariate") + geom_tile() + theme_pubr(x.text.angle = 90)
  return(p)
}
pdf("plots/covariate_correlation.pdf")
plot_heatmap(pc_tib)
dev.off()

# Plotting QQ plots of residual expression correlated with genotype PCs
pval_df = read_tsv("correlation.expression_regressed.global_ancestry.local_ancestry.race_Afr.race_Eur.seq_center.exam.genotype_PC1.txt")
QQ_plots = list()
i = 1
for (col in colnames(pval_df %>% select(-gene, -contains("expression")))) {
  tib = pval_df %>% transmute(pval = (!!as.name(col))) %>% arrange(pval) %>% mutate(expected = ppoints(pval))
  p = ggplot(tib, aes(-log10(expected), -log10(pval))) + geom_point() + 
    geom_abline(aes(intercept=0, slope=1)) + theme_pubr() + labs(title = col)
  QQ_plots[[i]] = p
  i = i + 1
}
pdf("plots/residual_expression.covariate_QQ_grid.pdf", width=15, height=10)
ggarrange(plotlist = QQ_plots, ncol=4, nrow=4)
dev.off()
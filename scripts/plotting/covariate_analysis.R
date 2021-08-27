library(tidyverse)
library(ggpubr)

# Set working directory to local directory containing paper figures
setwd("Google Drive/stanford/admixed PGS/figures/")
dir.create("plots", showWarnings = FALSE)

# Read covariate file
pc_tib = read_tsv("~/sherlock/oak/admixedPGS/eQTLs/results/covariate_analysis/v2/mesa.sample_covariates_deidentified.txt")

# Plot first two genotype PCs and color by race
tib = bind_rows(pc_tib %>% filter(race == 1) %>% mutate(race = "African-American"),
                pc_tib %>% filter(race == 0) %>% mutate(race = "European-American"))
ggplot(tib, aes(genotype_PC0, genotype_PC1, color = race)) + geom_point() + 
  theme_pubr() + theme(text = element_text(size=20)) + xlab("Genotype PC1") +
  scale_color_manual(values=c("#e48671", "#455A8E")) +
  ylab("Genotype PC2") +
  theme(axis.text.x = element_text(size=25), 
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5),
        legend.title = element_text(size=25),
        legend.text = element_text(size=25),
        legend.position = "bottom") +
  ggtitle("a)") + xlim(-160, 310) +
  guides(color = guide_legend(override.aes = list(size = 3) ) )
ggsave("plots/genotype_PCA.race.pdf", height=10, width=10, units="in")

# Plot correlation between genotype PC 1 and global African ancestry fraction
tib = pc_tib %>% filter(global_ancestry > 0)
ggplot(tib, aes(genotype_PC0, global_ancestry)) + geom_point() + 
  theme_pubr() + xlab("Genotype PC1") + ylab("African global ancestry fraction") + 
  stat_cor(digits=3, size=10) +
  theme(axis.text.x = element_text(size=25), 
        axis.text.y = element_text(size=25), 
        axis.title = element_text(size=30),
        plot.title = element_text(size=30),
        axis.line = element_line(size=0.5),
        axis.ticks = element_line(size=0.5)) +
  ggtitle("b)")
ggsave("plots/genotypePC0_vs_ganc.pdf", height=10, width=10, units="in")

# Plot QQ plots of expression correlated with covariates
pval_df = read_tsv("correlation.expression.txt")
QQ_plots = list()
i = 1
for (col in colnames(pval_df %>% select(-gene, -contains("expression")))) {
  tib = pval_df %>% transmute(pval = (!!as.name(col))) %>% arrange(pval) %>% mutate(expected = ppoints(pval))
  if (grepl("PC", col)) { # Convert 0-indexed PCs to 1-indexed PCs
    col = paste0(substr(col, 1, nchar(col) - 1), as.numeric(substr(col, nchar(col), nchar(col))) + 1)
  }
  p = ggplot(tib, aes(-log10(expected), -log10(pval))) + geom_point() + 
    geom_abline(aes(intercept=0, slope=1)) + theme_pubr() + labs(title = col) +
    theme(plot.title = element_text(size=20))
  QQ_plots[[i]] = p
  i = i + 1
}
ggarrange(plotlist = QQ_plots, ncol=4, nrow=4)
ggsave("plots/expression_covariate_QQ_grid.png", width=15, height=10, units="in")

# Plot expression PCs and color by batch
pdf("plots/expression_PCA.exam.pdf")
tib = bind_rows(pc_tib %>% filter(exam == 1) %>% mutate(exam = "Exam 1"),
                pc_tib %>% filter(exam == 0) %>% mutate(exam = "Exam 5"))
ggplot(tib, aes(expression_PC0, expression_PC1, color = exam)) + geom_point() + 
  theme_pubr() + theme(text = element_text(size=20)) + xlab("Expression PC1") +
  ylab("Expression PC2") + 
  scale_color_manual(values=c("#A396B1", "#F5A76C"))
dev.off()
pdf("plots/expression_PCA.seq_center.pdf")
tib = bind_rows(pc_tib %>% filter(seq_center == 1) %>% mutate(seq_center = "UW"),
                pc_tib %>% filter(seq_center == 0) %>% mutate(seq_center = "Broad"))
ggplot(tib, aes(expression_PC1, expression_PC2, color = seq_center)) + geom_point() + 
  theme_pubr() + theme(text = element_text(size=20)) + xlab("Expression PC2") +
  ylab("Expression PC3") + 
  scale_color_manual(values=c("#A396B1", "#F5A76C"))
dev.off()

# Plot correlation heatmap for covariates
plot_heatmap = function(tib) {
  cor_tib = as_tibble(tib %>% cor)
  long_cor_tib = cor_tib %>% mutate(covariate = colnames(cor_tib)) %>% pivot_longer(-covariate) %>% 
    mutate(value = abs(value))
  p = ggplot(long_cor_tib, aes(name, covariate, fill=value)) + scale_fill_distiller() + 
    xlab("covariate") + geom_tile() + theme_pubr(x.text.angle = 90)
  return(p)
}
pdf("plots/covariate_correlation.pdf")
plot_heatmap(pc_tib)
dev.off()

# Plot scree plot for genotype and expression PCs
plot_scree = function(filepath, n_eig) {
  eig = read_tsv(filepath, col_names =c("Eigenvalue"))
  eig = eig %>% mutate(Component = seq(1, nrow(eig)))
  p = ggplot(eig %>% head(n_eig), aes(Component, Eigenvalue)) + 
    geom_point() + theme_pubr() + theme(text = element_text(size=20))
  return(p)
}
pdf("plots/genotype_PC.top30_scree_plot.pdf")
plot_scree("~/sherlock/oak/admixedPGS/eQTLs/results/covariate_analysis/v2/mesa.prune.maf05.vcf.eigenvalues.txt", 30)
dev.off()
pdf("plots/expression_PC.top30_scree_plot.pdf")
plot_scree("~/sherlock/oak/admixedPGS/eQTLs/results/covariate_analysis/v2/mesa.expression.bed.eigenvalues.txt", 30)
dev.off()

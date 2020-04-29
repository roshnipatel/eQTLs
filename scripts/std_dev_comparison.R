library(tidyverse)
library(matrixStats)
library(optparse)

# Parse arguments
option_list = list(make_option("--eur", type="character", default=NULL,
                        help="European expression file name", metavar="character"),
                   make_option("--afr", type="character", default=NULL,
                        help="African expression file name", metavar="character"),
                   make_option("--samp_prefix", type="character", default=NULL,
                        help="filepath prefix for sample input files", metavar="character"),
                   make_option("--out", type="character", default=NULL,
                        help="output filepath", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

afr = read_tsv(opt$afr)
eur = read_tsv(opt$eur)
all_ind = afr %>% inner_join(eur)

# Within each population, compute standard deviation in gene reads for all genes
afr_sd = afr %>% select(-'#chr', -start, -end, -gene_id) %>% transmute(sd = rowSds(as.matrix(.), na.rm=TRUE)) %>% bind_cols(afr %>% select(gene_id))
eur_sd = eur %>% select(-'#chr', -start, -end, -gene_id) %>% transmute(sd = rowSds(as.matrix(.), na.rm=TRUE)) %>% bind_cols(eur %>% select(gene_id))
exp_sd = afr_sd %>% inner_join(eur_sd, by="gene_id", suffix=c("_afr", "_eur"))

# Define function to find expression SD in pop, a subset of AA individuals, for gene
pop_subset_sd = function(pop, gene) {
    gene_file = paste0(args$samp_prefix, pop, "/", gene, ".txt")
    if (!file.exists(gene_file)) {
        return(NA)
    }
    ind = read_tsv(gene_file, col_names=c("NWDID"))
    row_sd = all_ind %>% filter(gene_id == gene) %>% select(ind$NWDID) %>% sd(na.rm=TRUE)
    return(row_sd)
}
pop_subset_sd_vect = Vectorize(pop_subset_sd)

exp_sd = exp_sd %>% mutate(Afr_sd = pop_subset_sd_vect("Afr", gene_id), het_sd = pop_subset_sd_vect("het", gene_id), Eur_sd = pop_subset_sd_vect("Eur", gene_id))
exp_sd %>% write.table(args$out, row.names=FALSE, sep='\t')

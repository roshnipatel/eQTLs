library(tidyverse)
library(optparse)
library(qvalue)

# Parse arguments
option_list = list(make_option("--asc", type="character", default=NULL,
                        help="merged ascertainment file name", metavar="character"),
                   make_option("--afr", type="character", default=NULL,
                        help="merged afr estimation file name", metavar="character"),
                   make_option("--eur", type="character", default=NULL,
                        help="merged eur estimation file name", metavar="character"),
                   make_option("--het", type="double", default=NULL,
                        help="merged het estimation file name", metavar="character"),
                   make_option("--out", type="character", default=NULL,
                        help="output file name prefix", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

asc = read_tsv(opt$asc)

# Apply Storey-Tibshirani FDR correction  to ascertainment dataset to identify significant hits
asc = asc %>% mutate(fdr = qvalue(perm_pval)$qvalues)
hits = asc %>% filter(fdr <= .10)
hits %>% write.table(paste0(opt$out, "hits_ascertainment_Eur.txt"), sep='\t', row.names=FALSE)

# Extract effect sizes (+ other info) for significant hits from estimation datasets
# Note that estimation datasets must be read and filtered sequentially otherwise you'll run out of memory
eur = read_tsv(opt$eur)
eur = eur %>% inner_join(hits %>% select(gene, variant))
eur %>% write.table(paste0(opt$out, "hits_estimation_Eur.txt", sep='\t', row.names=FALSE)

afr = read_tsv(opt$afr)
afr = afr %>% inner_join(hits %>% select(gene, variant))
afr %>% write.table(paste0(opt$out, "hits_estimation_Afr.txt", sep='\t', row.names=FALSE)

het = read_tsv(opt$het)
het = het %>% inner_join(hits %>% select(gene, variant))
het %>% write.table(paste0(opt$out, "hits_estimation_het.txt", sep='\t', row.names=FALSE)

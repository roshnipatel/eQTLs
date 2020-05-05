library(tidyverse)
library(optparse)
library(qvalue)

# Parse arguments
option_list = list(make_option("--asc", type="character", default=NULL,
                        help="merged ascertainment file name", metavar="character"),
                   make_option("--hit", type="character", default=NULL,
                        help="ascertainment hits file name", metavar="character"),
                   make_option("--afr", type="character", default=NULL,
                        help="merged afr estimation file name", metavar="character"),
                   make_option("--eur", type="character", default=NULL,
                        help="merged eur estimation file name", metavar="character"),
                   make_option("--het", type="character", default=NULL,
                        help="merged het estimation file name", metavar="character"),
                   make_option("--val", type="character", default=NULL,
                        help="merged eur validation file name", metavar="character"),
                   make_option("--val_afr", type="character", default=NULL,
                        help="merged afr validation file name", metavar="character"),
                   make_option("--out", type="character", default=NULL,
                        help="output file name prefix", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

USING_PERM = FALSE

if (is.null(opt$asc) & !is.null(opt$hit)) {
    hits = read_tsv(opt$hit)
} else if (is.null(opt$hit) & !is.null(opt$asc)) {
    asc = read_tsv(opt$asc)
    
    if (USING_PERM) {
        # Apply Storey-Tibshirani FDR correction  to ascertainment dataset to identify significant hits
        asc = asc %>% mutate(fdr = qvalue(perm_pval)$qvalues)
        hits = asc %>% filter(fdr <= .10)
        hits %>% write.table(paste0(opt$out, "hits_ascertainment_Eur.txt"), sep='\t', row.names=FALSE)
    } else {
        asc = asc %>% mutate(fdr = qvalue(pval)$qvalues)
        hits = asc %>% filter(fdr < .10) %>% group_by(Gene) %>% filter(fdr == min(fdr), pval == min(pval)) %>% sample_n(1) %>% ungroup
        hits %>% write.table(paste0(opt$out, "hits_ascertainment_Eur.txt"), sep='\t', row.names=FALSE)
    }
}

# Extract effect sizes (+ other info) for significant hits from estimation datasets
# Note that estimation datasets must be read and filtered sequentially otherwise you'll run out of memory
eur = read_tsv(opt$eur)
eur = eur %>% inner_join(hits %>% select(Gene, ID))
eur %>% write.table(paste0(opt$out, "hits_estimation_Eur.txt"), sep='\t', row.names=FALSE)

afr = read_tsv(opt$afr)
afr = afr %>% inner_join(hits %>% select(Gene, ID))
afr %>% write.table(paste0(opt$out, "hits_estimation_Afr.txt"), sep='\t', row.names=FALSE)

het = read_tsv(opt$het)
het = het %>% inner_join(hits %>% select(Gene, ID))
het %>% write.table(paste0(opt$out, "hits_estimation_het.txt"), sep='\t', row.names=FALSE)

val = read_tsv(opt$val)
val = val %>% inner_join(hits %>% select(Gene, ID))
val %>% write.table(paste0(opt$out, "hits_validation_Eur.txt"), sep='\t', row.names=FALSE)

val_afr = read_tsv(opt$val_afr)
val_afr = val_afr %>% inner_join(hits %>% select(Gene, ID))
val_afr %>% write.table(paste0(opt$out, "hits_validation_Afr.txt"), sep='\t', row.names=FALSE)

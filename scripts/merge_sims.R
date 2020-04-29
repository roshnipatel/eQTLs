library(tidyverse)
library(optparse)

# Parse arguments
option_list = list(make_option("--dir", type="character", default=NULL,
                        help="directory with simulation files", metavar="character"),
                   make_option("--out", type="character", default=NULL,
                        help="output file name", metavar="character"))
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

sim_files = paste0(opt$dir, list.files(opt$dir))

sim_res = tibble()
for (f in sim_files) {
    tmp = read_tsv(f)
    sim_res = sim_res %>% bind_rows(tmp)
}
sim_res %>% write.table(opt$out, sep='\t')

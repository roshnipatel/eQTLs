library(tidyverse)

reads = read_tsv("data/TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_reads.gct.gz", skip=2)
metadata = read_csv("data/MESA_sample_info.csv")
samples = read_tsv("data/sample_participant_lookup.txt")

AA = samples %>% left_join(metadata %>% select(NWDID, race1c) %>% unique) %>% filter(race1c = 3)
EA = samples %>% left_join(metadata %>% select(NWDID, race1c) %>% unique) %>% filter(race1c = 1)
AAreads = reads %>% select(AA$TOR_ID)
EAreads = reads %>% select(EA$TOR_ID)
AAsd = AAreads %>% transmute(AA_sd = rowSds(as.matrix(.[AAreads %>% colnames])))
EAsd = EAreads %>% transmute(EA_sd = rowSds(as.matrix(.[EAreads %>% colnames])))

expsd = bind_cols(reads %>% select(Name, Description), AAsd, EAsd)
expsd %>% write.table("data/expression_std_dev.txt", row.names=FALSE, sep='\t')

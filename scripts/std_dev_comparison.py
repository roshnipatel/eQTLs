import pandas as pd
import scipy.stats
import numpy as np

def compare_var(df1, df2, sd=False):
    comb = pd.merge(df1.var(axis=1).rename("one"), df2.var(axis=1).rename("two"), left_index=True, right_index=True)
    if sd:
        res = scipy.stats.ttest_rel(np.sqrt(comb["one"]), np.sqrt(comb["two"]))
    else:
        res = scipy.stats.ttest_rel(comb["one"], comb["two"])
    return(res)

if __name__ == "__main__":
    afr_props = pd.read_csv("data/mesa.Afr.expression.bed.gz", sep='\t').set_index("gene_id").drop(["#chr", "start", "end"], axis=1)
    eur_props = pd.read_csv("data/mesa.Eur.expression.bed.gz", sep='\t').set_index("gene_id").drop(["#chr", "start", "end"], axis=1)
    props = pd.merge(afr_props.var(axis=1).rename("afr_props_var"), eur_props.var(axis=1).rename("eur_props_var"), left_index=True, right_index=True)
    
    reads = pd.read_csv("data/TOPMed_MESA_RNAseq_Pilot_expression_data/TOPMed_MESA_RNAseq_Pilot_RNASeQCv1.1.9.gene_reads.gct.gz", sep='\t', skiprows=2).set_index("Name").drop("Description", axis=1)
    
    afr_samp = pd.read_csv("data/select_samples_Afr.txt", sep='\t')
    eur_samp = pd.read_csv("data/select_samples_Eur.txt", sep='\t')
    
    afr_reads = reads[afr_samp.tor_id].loc[afr_props.index]
    eur_reads = reads[eur_samp.tor_id].loc[eur_props.index]
    reads = pd.merge(afr_reads.var(axis=1).rename("afr_reads_var"), eur_reads.var(axis=1).rename("eur_reads_var"), left_index=True, right_index=True)
    
    all_var = pd.merge(props, reads, left_index=True, right_index=True)
    all_var.to_csv("data/variances.txt", sep='\t')

    # props_res = compare_var(afr_props, eur_props)
    # print(props_res)
    # reads_res = compare_var(afr_reads, eur_reads)
    # print(reads_res)
    # props_res = compare_var(afr_props, eur_props, True)
    # print(props_res)
    # reads_res = compare_var(afr_reads, eur_reads, True)
    # print(reads_res)

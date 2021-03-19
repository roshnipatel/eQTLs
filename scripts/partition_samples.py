import argparse
import pandas as pd
import numpy as np

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--intersect')
    parser.add_argument('--afr_samples')
    parser.add_argument('--eur_samples')
    parser.add_argument('--genes')
    parser.add_argument('--ascertainment')
    parser.add_argument('--validation')
    parser.add_argument('--ascertainment_pop', default="Eur", help="either 'Afr' or 'Eur'")
    parser.add_argument('--afr_validation_run', action='store_true')
    parser.add_argument('--out_dir')
    args = parser.parse_args()
    
    CIS_WINDOW = 200000
    
    afr_samples = pd.read_csv(args.afr_samples, delimiter='\t')
    eur_samples = pd.read_csv(args.eur_samples, delimiter='\t')
    genes = pd.read_csv(args.genes, delimiter='\t', 
                        names=["chrom", "start", "stop", "gene_id"])["gene_id"]
    
    # Filter tracts data for individuals and genes in our sample dataset.
    tracts = pd.read_csv(args.intersect, delimiter='\t',
                        names=["chrom", "anc_start", "anc_stop", "info", "gene_chrom", 
                               "gene_start", "gene_stop", "gene_id", "overlap"],
                        dtype={"chrom": int, "anc_start": int, "anc_stop": int, 
                               "info": str, "gene_chrom": str, "gene_start": int, 
                               "gene_stop": int, "gene_id": str, "overlap": int})
    tracts = pd.merge(genes, tracts, how='inner')
    tracts["nwd_id"] = tracts.apply(lambda row: row.info.split('_')[0], axis=1)
    tracts["anc"] = tracts.apply(lambda row: row.info.split('_')[2], axis=1)
    tracts = pd.merge(afr_samples, tracts, how='left')[["gene_id", "nwd_id", 
                                                        "anc", "overlap"]]
    
    # Identify tracts of African ancestry
    Afr_tracts = tracts[(tracts["overlap"] == CIS_WINDOW) & (tracts["anc"] == 'YRI')]
    Afr_tracts = Afr_tracts.drop(['anc'], axis=1)
    Afr_tracts = Afr_tracts.groupby(["gene_id", "nwd_id"]).size()
    
    # Identify tracts of European ancestry
    Eur_tracts = tracts[(tracts["overlap"] == CIS_WINDOW) & (tracts["anc"] == 'CEU')]
    Eur_tracts = Eur_tracts.drop(["anc"], axis=1)
    Eur_tracts = Eur_tracts.groupby(["gene_id", "nwd_id"]).size()
    
    # Identify Afr-Am individuals that have African ancestry overlapping gene's
    # cis window for both chromosomes.
    hom_Afr_tracts = Afr_tracts[Afr_tracts > 1]
    hom_Afr_tracts = hom_Afr_tracts.reset_index(name='counts')
    hom_Afr_tracts = hom_Afr_tracts.drop(['counts'], axis=1)

    # Initialize matrices to store ascertainment and validation information. 
    # Note that this validation matrix doesn't store the IDs of individuals in 
    # the directory 'reestimation_validation', which is used for results in the 
    # QTL_calling directory. Instead, this matrix stores a random sample of 20% 
    # of individuals to be used for variance partitioning analyses.
    all_ind = pd.concat([eur_samples.nwd_id, afr_samples.nwd_id])
    all_genes = gene_counts.keys()
    def create_matrix(genes, ind):
        mat = pd.DataFrame(np.zeros((len(genes), len(ind))))
        mat.columns = ind
        mat.index = genes
        return(mat)
    ascertainment_matrix = create_matrix(all_genes, all_ind)
    validation_matrix = create_matrix(all_genes, all_ind)
    
    # Identify Afr-Am individuals that have one chromosome with African ancestry
    # overlapping gene's cis window and one chromosome with European ancestry 
    # overlapping gene's cis window.
    het_Afr_tracts = Afr_tracts[Afr_tracts == 1]
    het_Afr_tracts = het_Afr_tracts.reset_index(name='counts')
    het_Eur_tracts = Eur_tracts[Eur_tracts == 1]
    het_Eur_tracts = het_Eur_tracts.reset_index(name='counts')
    anc_het_tracts = pd.merge(het_Afr_tracts, het_Eur_tracts, how='inner')
    
    # Write ancestry-heterozygous sample IDs to file
    for gene, df in anc_het_tracts.groupby("gene_id"):
        filepath = args.out_dir + "/reestimation_primary/het/" + gene + ".txt"
        df.to_csv(filepath, header=False, index=False, columns=["nwd_id"])
        
    # Optionally randomly samply 50 African-Americans to serve as validation set
    # for sanity-checking QTL_calling results. (Different from the validation 
    # matrix created for variance partitioning analyses!)
    if args.afr_validation_run:
        for gene, df in hom_Afr_tracts.groupby("gene_id"):
            filepath = args.out_dir + "/reestimation_validation/Afr/" + gene + ".txt"
            df.sample(n = 50).to_csv(filepath, header=False, 
                                     index=False, columns=["nwd_id"])
   
    # Partition ascertainment/estimation individuals depending on desired 
    # ascertainment population.
    if args.ascertainment_pop == "Eur":    
        # Write Afr-Am sample IDs to file and store sample sizes for each gene
        gene_counts = {}
        for gene, df in hom_Afr_tracts.groupby("gene_id"):
            gene_counts[gene] = df.shape[0]
            filepath = args.out_dir + "/reestimation_primary/Afr/" + gene + ".txt"
            df.to_csv(filepath, header=False, index=False, columns=["nwd_id"])  
        
        Eur_IDs = eur_samples["nwd_id"]
        for gene, count in gene_counts.items():
            # Randomly sample same number of European-Americans to estimate
            # SNP effect sizes.
            est = Eur_IDs.sample(n = count)
            est_filepath = args.out_dir + "/reestimation_primary/Eur/" + gene + ".txt"
            est.to_csv(est_filepath, header=False, index=False)

            # Randomly sample 50 European-Americans to be in estimation validation
            # set (deprecated part of analysis pipeline).
            non_est = Eur_IDs.drop(est.index)
            val = non_est.sample(n = 50)
            val_filepath = args.out_dir + "/reestimation_validation/Eur/" + gene + ".txt"
            val.to_csv(val_filepath, header=False, index=False)

            # Remaining individuals comprise ascertainment set.
            asc = non_est.drop(val.index)
            ascertainment_matrix.loc[gene, asc] = 1
            asc_filepath = args.out_dir + "/ascertainment/Eur/" + gene + ".txt"
            asc.to_csv(asc_filepath, header=False, index=False)

            # Randomly sample 20% of non-ascertainment individuals to be in
            # validation set for variance partitioning analyses.
            non_asc = Eur_IDs.drop(asc.index)
            n_val_anova = int(len(non_asc) * .2)
            val_anova = non_asc.sample(n = n_val_anova)
            validation_matrix.loc[gene, val_anova] = 1
        ascertainment_matrix.to_csv(args.ascertainment, sep='\t', header=True, index=True)
        
        # Randomly sample 20% of African-Americans to be in validation set for
        # variance partitioning analyses.
        Afr_IDs = afr_samples["nwd_id"]
        n_Afr_val = int(len(Afr_IDs) * .2)
        Afr_val = Afr_IDs.sample(n = n_Afr_val)
        validation_matrix.loc[all_genes, Afr_val] = 1
        validation_matrix.to_csv(args.validation, sep='\t', header=True, index=True)
    elif args.ascertainment_pop == "Afr":
        gene_counts = {}
        for gene, df in hom_Afr_tracts.groupby("gene_id"):
            # Randomly sample 2/3 of African-Americans with 2 copies of African
            # ancestry to be ascertainment set.
            IDs = df["nwd_id"]
            n_asc = int(IDs.shape[0] * 2 / 3)
            asc = IDs.sample(n = n_asc)
            ascertainment_matrix.loc[gene, asc] = 1
            asc_filepath = args.out_dir + "/ascertainment/Afr/" + gene + ".txt"
            asc.to_csv(asc_filepath, header=False, index=False)

            # Remaining African-Americans used for (re)estimation of SNP effect 
            # size.
            est = IDs.drop(asc.index)
            gene_counts[gene] = est.shape[0]
            est_filepath = args.out_dir + "/reestimation_primary/Afr/" + gene + ".txt"
            est.to_csv(est_filepath, header=False, index=False)

            # Randomly sample 20% of non-ascertainment individuals to be in
            # validation set for variance partitioning analyses.
            n_val_anova = int(len(est) * .2)
            val_anova = est.sample(n = n_val_anova)
            validation_matrix.loc[gene, val_anova] = 1
        ascertainment_matrix.to_csv(args.ascertainment, sep='\t', header=True, index=True)
        
        est = eur_samples["nwd_id"]
        for gene, count in gene_counts.items():
            # Randomly sample same number of European-Americans to estimate
            # SNP effect sizes in VALIDATION SET (this is different from what
            # we do for European-based ascertainment).
            val = est.sample(n = count)
            val_filepath = args.out_dir + "/reestimation_validation/Eur/" + gene + ".txt"
            val.to_csv(val_filepath, header=False, index=False)

            # Use ALL European-Americans to comprise estimation set.
            est_filepath = args.out_dir + "/reestimation_primary/Eur/" + gene + ".txt"
            est.to_csv(est_filepath, header=False, index=False)

        # Randomly sample 20% of European-Americans to be in validation set for
        # variance partitioning analyses.
        n_Eur_val = int(len(est) * .2)
        Eur_val = est.sample(n = n_Eur_val)
        validation_matrix.loc[all_genes, Eur_val] = 1
        validation_matrix.to_csv(args.validation, sep='\t', header=True, index=True)


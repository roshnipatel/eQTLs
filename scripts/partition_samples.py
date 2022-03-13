"""
This script creates 6 different types of sample input files.
  -> ascertainment/Eur stores random sample of 232 EA for each gene.
  -> reestimation/primary/EA stores remaining 250 EA for each gene.
  -> reestimation/primary/AA stores random sample of 250 AA for each gene.
  -> reestimation/primary/Afr stores AA with 2 Afr haplo for each gene (max 250).
  -> reestimation/primary/Eur stores random sample of same number of non-asc EA.
  -> reestimation/primary/het stores AA with 1 Afr haplo for each gene.
args.ascertainment stores a matrix of ascertainment (EA) individuals.
args.validation stores a random sample of 20% of AA and EA to be used as a 
validation set when computing variance explained by different models.
"""

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
    parser.add_argument('--out_dir')
    args = parser.parse_args()
    
    CIS_WINDOW = 200000
    MAX_AA = 250
    
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
    
    # Initialize matrices to store ascertainment and validation information. 
    # Note that this validation matrix doesn't store the IDs of individuals in 
    # the directory 'reestimation_validation', which is used for results in the 
    # QTL_calling directory. Instead, this matrix stores a random sample of 20% 
    # of individuals to be used for variance partitioning analyses.
    all_ind = pd.concat([eur_samples.nwd_id, afr_samples.nwd_id])
    def create_matrix(genes, ind):
        mat = pd.DataFrame(np.zeros((len(genes), len(ind))))
        mat.columns = ind
        mat.index = genes
        return(mat)
    ascertainment_matrix = create_matrix(genes, all_ind)
    validation_matrix = create_matrix(genes, all_ind)
    
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

    # # Identify Afr-Am individuals that have one chromosome with African ancestry
    # # overlapping gene's cis window and one chromosome with European ancestry 
    # # overlapping gene's cis window. Deprecated.
    # het_Afr_tracts = Afr_tracts[Afr_tracts == 1]
    # het_Afr_tracts = het_Afr_tracts.reset_index(name='counts')
    # het_Eur_tracts = Eur_tracts[Eur_tracts == 1]
    # het_Eur_tracts = het_Eur_tracts.reset_index(name='counts')
    # anc_het_tracts = pd.merge(het_Afr_tracts, het_Eur_tracts, how='inner')
    # for gene, df in anc_het_tracts.groupby("gene_id"):
    #     filepath = args.out_dir + "/reestimation_primary/het/" + gene + ".txt"
    #     df.to_csv(filepath, header=False, index=False, columns=["nwd_id"])
        
    # Partition ascertainment/estimation individuals depending on desired 
    # ascertainment population.
    if args.ascertainment_pop == "Eur":    
        n_validation_AA = int(len(afr_samples.nwd_id) * .2)
        n_validation_EA = int(MAX_AA * .2)

        for gene in genes:
            if gene == "gene_id":
                continue
            # African-American estimation set (indep. of local ancestry)
            AA_est = afr_samples.nwd_id.sample(n = MAX_AA)
            AA_filepath = args.out_dir + "/reestimation_primary/AA/" + gene + ".txt"
            AA_est.to_csv(AA_filepath, header=False, index=False)

            # African-American validation set
            AA_val = afr_samples.nwd_id.sample(n = n_validation_AA)
            validation_matrix.loc[gene, AA_val] = 1

            # African estimation set (2 haplotypes African ancestry)
            Afr_est = hom_Afr_tracts[hom_Afr_tracts.gene_id == gene].nwd_id
            n_Afr_est = len(Afr_est)
            if n_Afr_est > MAX_AA:
                Afr_est = Afr_est.sample(n = MAX_AA)
                n_Afr_est = len(Afr_est)
            Afr_filepath = args.out_dir + "/reestimation_primary/Afr/" + gene + ".txt"
            Afr_est.to_csv(Afr_filepath, header=False, index=False)

            # European-American estimation set (n = AA estimation set)
            EA_est = eur_samples.nwd_id.sample(n = MAX_AA)
            EA_filepath = args.out_dir + "/reestimation_primary/EA/" + gene + ".txt"
            EA_est.to_csv(EA_filepath, header=False, index=False)

            # European-American validation set
            EA_val = EA_est.sample(n = n_validation_EA)
            validation_matrix.loc[gene, EA_val] = 1

            # European estimation set (n = Afr estimation set)
            Eur_est = EA_est.sample(n = n_Afr_est)
            Eur_filepath = args.out_dir + "/reestimation_primary/Eur/" + gene + ".txt"
            Eur_est.to_csv(Eur_filepath, header=False, index=False)

            # European-American ascertainment set
            asc = eur_samples.nwd_id.drop(EA_est.index)
            ascertainment_matrix.loc[gene, asc] = 1
            asc_filepath = args.out_dir + "/ascertainment/Eur/" + gene + ".txt"
            asc.to_csv(asc_filepath, header=False, index=False)
    
        ascertainment_matrix.to_csv(args.ascertainment, sep='\t', header=True, index=True)
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


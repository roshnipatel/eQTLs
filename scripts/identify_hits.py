import argparse
import pandas as pd

def FDR_threshold(tag_SNPs, colname, fdr):
    """Apply FDR threshold to SNPs based on Benjamini-Hochberg procedure."""
    n_SNPs = tag_SNPs.shape[0]
    def check_significance(row):
        pval_threshold = (row.name + 1) / n_SNPs * fdr
        return(not (row[colname] > pval_threshold))
    ordered_SNPs = tag_SNPs.sort_values(by=colname, ignore_index=True)
    ordered_SNPs["significant"] = ordered_SNPs.apply(check_significance, axis=1)
    sig_SNPs = ordered_SNPs[ordered_SNPs.significant == True]
    return(sig_SNPs)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ascertainment', action='store_true')
    parser.add_argument('--reestimation', action='store_true')
    parser.add_argument('--merged')
    parser.add_argument('--tss')
    parser.add_argument('--hits', help="existing file of hits to select (used in reestimation mode)")
    parser.add_argument('--fdr', help="FDR threshold (used in ascertainment mode)")
    parser.add_argument('--out')
    args = parser.parse_args()

    if args.ascertainment: 
        # Identify hits in ascertainment dataset based on user-specified FDR.
        # For each gene, select SNP with smallest p-value if multiple pass the 
        # FDR threshold.
        fdr = float(args.fdr)
        all_SNPs = pd.read_csv(args.merged, sep='\t')
        gene_tss = pd.read_csv(args.tss, sep='\t', usecols=["start", "gene_id"])
        all_SNPs = pd.merge(all_SNPs, gene_tss, left_on="gene", right_on="gene_id")
        all_SNPs["pos"] = all_SNPs.apply(lambda row: int(row.ID.split('_')[1]), axis=1)
        all_SNPs["dist"] = all_SNPs.apply(lambda row: abs(row.pos - row.start), axis=1)
        all_SNPs["pval_round"] = all_SNPs.apply(lambda row: round(row.pval, 4), axis=1)

        # First select SNPs with smallest p-value for each gene
        min_pvals = all_SNPs.groupby("gene").pval_round.min().reset_index()
        min_pval_SNPs = pd.merge(all_SNPs, min_pvals)

        # Then choose nearest SNP to break ties
        min_pval_dist_SNPs = min_pval_SNPs.loc[min_pval_SNPs.groupby("gene").dist.idxmin()]
        sig_SNPs = FDR_threshold(min_pval_dist_SNPs, "pval", fdr)
        sig_SNPs.to_csv(args.out, index=False, sep='\t')
    elif args.reestimation: 
        # Extract hits from reestimation dataset based on existing ascertainment hits
        all_SNPs = pd.read_csv(args.merged, sep='\t')
        hits = pd.read_csv(args.hits, sep='\t')
        sig_SNPs = pd.merge(all_SNPs, hits[["gene", "ID"]])
        sig_SNPs.to_csv(args.out, index=False, sep='\t')

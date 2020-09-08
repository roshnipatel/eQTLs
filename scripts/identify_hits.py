import argparse
import pandas as pd

def FDR_threshold(tag_SNPs, colname, fdr):
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
    parser.add_argument('--merged', help="merged file path to extract hits from")
    parser.add_argument('--hits', help="existing file of hits to select (used in reestimation mode)")
    parser.add_argument('--out', help="output file path")
    parser.add_argument('--fdr', help="FDR threshold (used in ascertainment mode)")
    args = parser.parse_args()

    if args.ascertainment: # Identify hits in ascertainment dataset based on user-specified FDR
        fdr = float(args.fdr)
        all_SNPs = pd.read_csv(args.merged, sep='\t')
        sig_SNPs = FDR_threshold(all_SNPs, "pval", fdr)
        sig_SNPs = sig_SNPs.loc[sig_SNPs.groupby("gene")["pval"].idxmin()] # Select max one SNP per gene
        sig_SNPs.to_csv(args.out, index=False, sep='\t')
    elif args.reestimation: # Identify hits in reestimation dataset based on existing ascertainment hits
        all_SNPs = pd.read_csv(args.merged, sep='\t')
        hits = pd.read_csv(args.hits, sep='\t')
        sig_SNPs = pd.merge(all_SNPs, hits[["gene", "ID"]])
        sig_SNPs.to_csv(args.out, index=False, sep='\t')

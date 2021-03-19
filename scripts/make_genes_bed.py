import pandas as pd
import argparse

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--genes')
    parser.add_argument('--window')
    parser.add_argument('--chrom_lengths')
    parser.add_argument('--out')
    args = parser.parse_args()

    # Specify radius of cis window around TSS
    window = int(args.window)
    
    # Load list of genes that pass QC steps in normalization step
    genes = pd.read_csv(args.genes, delimiter='\t', skiprows=1,
                        names=["chrom", "pos_before_tss", "tss", "gene_id"])
    genes["chrom"] = genes.apply(lambda row: row.chrom[3:], axis=1)
    
    # Load list of chromosomes and their lengths
    chrom_lengths = pd.read_csv(args.chrom_lengths, delimiter='\t', 
                                names=["chrom", "length"])
    
    # Define start and end of cis window for each gene
    genes["cis_start"] = genes.apply(lambda row: max(row.tss - window, 0), axis=1)
    genes["cis_stop"] = genes.apply(lambda row: min(row.tss + window, 
        chrom_lengths[chrom_lengths["chrom"] == row.chrom]["length"].iloc[0]), axis=1)
    
    genes.to_csv(args.out, columns=["chrom", "cis_start", "cis_stop", "gene_id"], 
                 sep='\t', index=False, header=False)

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--anno', help='gene annotation file')
parser.add_argument('--chrom_lengths')
parser.add_argument('--genes')
parser.add_argument('--out', help='output')
args = parser.parse_args()

offset = 100000

genes = pd.read_csv(args.genes, delimiter='\t', skiprows=1, names=["chrom", "pos_before_tss", "tss", "gene_id"])[["gene_id", "tss"]]
chrom_lengths = pd.read_csv(args.chrom_lengths, delimiter='\t', names=["chrom", "length"])

anno = pd.read_csv(args.anno, delimiter='\t', skiprows=6, names=["chrom", "transcript", "type", "gene_start", "gene_stop", "x1", "x2", "x3", "info"])[["chrom", "info", "type"]]
anno = anno[anno.type == 'gene'] # Ignore exons and transcripts
anno["gene_id"] = anno.apply(lambda row: row.info.split(';')[0].split('"')[1], axis=1) # Extract Ensembl gene ID from the semicolon-delimited info field
anno = pd.merge(genes, anno, how='left') # Subset for genes that we have expression data for
anno["chrom"] = anno.apply(lambda row: row.chrom[3:], axis=1)
anno["cis_start"] = anno.apply(lambda row: max(row.tss - offset, 0), axis=1)
anno["cis_stop"] = anno.apply(lambda row: min(row.tss + offset, chrom_lengths[chrom_lengths["chrom"] == row.chrom]["length"].iloc[0]), axis=1)

anno.to_csv(args.out, columns=["chrom", "cis_start", "cis_stop", "gene_id"], sep='\t', index=False, header=False)

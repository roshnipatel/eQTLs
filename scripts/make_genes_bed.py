import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--anno', help='gene annotation file')
parser.add_argument('--chrom_lengths')
parser.add_argument('--genes')
parser.add_argument('--out', help='output')
args = parser.parse_args()

offset = 1000000

genes = pd.read_csv(args.genes, delimiter='\t', names=["Chrom", "Pos_before_TSS", "TSS", "GeneID"])[["GeneID", "TSS"]]
chrom_lengths = pd.read_csv(args.chrom_lengths, delimiter='\t', names=["Chrom", "Length"])

anno = pd.read_csv(args.anno, delimiter='\t', skiprows=6, names=["Chrom", "Transcript", "Type", "GeneStart", "GeneStop", "x1", "x2", "x3", "Info"])[["Chrom", "Info", "Type"]]
anno = anno[anno.Type == 'gene'] # Ignore exons and transcripts
anno["GeneID"] = anno.apply(lambda row: row.Info.split(';')[0].split('"')[1], axis=1) # Extract Ensembl gene ID from the semicolon-delimited Info field
anno = pd.merge(genes, anno, how='left') # Subset for genes that we have expression data for
anno["Chrom"] = anno.apply(lambda row: row.Chrom[3:], axis=1)
anno["cisStart"] = anno.apply(lambda row: max(row.TSS - offset, 0), axis=1)
anno["cisStop"] = anno.apply(lambda row: min(row.TSS + offset, chrom_lengths[chrom_lengths["Chrom"] == row.Chrom]["Length"].iloc[0]), axis=1)

anno.to_csv(args.out, columns=["Chrom", "cisStart", "cisStop", "GeneID"], sep='\t', index=False, header=False)

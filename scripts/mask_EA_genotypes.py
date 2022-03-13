import pandas as pd
import numpy as np
import gzip
import argparse
from mask_AA_genotypes import parse_genotypes, drop_genotypes, write_vcf

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf')
    parser.add_argument('--out', nargs='+')
    args = parser.parse_args()

    vcf_header, genotypes, geno_ind = parse_genotypes(args.vcf)
    swapped_geno = genotypes.copy()
    # Randomly mask 1/3 of individuals
    for ind in geno_ind:
        x = np.random.uniform(0, 3)
        if x > 2:
            genotypes.loc[:,ind] = genotypes.loc[:,ind].apply(drop_genotypes)
        else:
            swapped_geno.loc[:,ind] = swapped_geno.loc[:,ind].apply(drop_genotypes)

    write_vcf(vcf_header, genotypes, args.out[0])
    write_vcf(vcf_header, swapped_geno, args.out[1])

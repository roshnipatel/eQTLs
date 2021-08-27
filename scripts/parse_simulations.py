import pandas as pd
import numpy as np
import argparse
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--delta_dir')
    parser.add_argument('--out', help='output file')
    args = parser.parse_args()

    delta_files = [os.path.join(args.delta_dir, f) for f in os.listdir(args.delta_dir)]
    sim_results = []
    for path in delta_files:
        with open(path, 'r') as f:
            estimated_delta = float(f.readlines()[-1])
        true_delta = float(path[-26:-23])
        sim_results.append([true_delta, estimated_delta])
   
    sim_results = pd.DataFrame(sim_results, columns=["true_delta", "estimated_delta"])
    sim_results.to_csv(args.out, index=False, sep='\t')

import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--delta_files', nargs='+')
    parser.add_argument('--simulated_delta')
    parser.add_argument('--out', help='output file')
    args = parser.parse_args()

    sim_delta = float(args.simulated_delta)

    delta_vals = []
    for path in args.delta_files:
        with open(path, 'r') as f:
            final_delta = float(f.readlines()[-1])
        delta_vals.append(final_delta)
    delta_vals = np.array(delta_vals)
    mean_delta = np.mean(delta_vals)
    bias_delta = mean_delta - sim_delta
    var_delta = np.var(delta_vals)
    
    with open(args.out, 'w') as f:
        f.write('Bias: {0}\n'.format(bias_delta))
        f.write('Variance: {0}\n'.format(var_delta))

import pandas as pd
import numpy as np
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--bootstrap_files', nargs='+')
    parser.add_argument('--out', help='prefix for output files')
    args = parser.parse_args()

    # Combine estimates of delta from all bootstrap iterations
    delta_vals = []
    for path in args.bootstrap_files:
        print(path)
        with open(path, 'r') as f:
            final_delta = float(f.readlines()[-1])
            print(final_delta)
        delta_vals.append(final_delta)
    delta_vals = np.array(delta_vals)
    lower = np.quantile(delta_vals, .025)
    upper = np.quantile(delta_vals, .975)
    
    # Write 95% confidence interval for delta to output
    with open(args.out + "summary.txt", 'w') as f:
        f.write('.025 quantile: {0}\n'.format(lower))
        f.write('.975 quantile: {0}\n'.format(upper))
    np.savetxt(args.out + "all.txt", delta_vals, delimiter='\n')

import argparse

import pandas as pd
import numpy as np


def parse_cl_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input', dest='input_file', required=True)
    parser.add_argument('-o', '--output', dest='output', required=True)
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_cl_args()

    dge_df = pd.read_table(args.input_file, index_col=0)
    dge_df['rnk'] = -np.log(dge_df['pvalue']) * dge_df['log2FoldChange']
    rnk_df = dge_df['rnk'].dropna().sort_values(ascending=False)
    rnk_df.to_csv(args.output, sep='\t', header=False)
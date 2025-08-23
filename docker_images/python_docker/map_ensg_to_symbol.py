import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input_matrix')
    parser.add_argument('-m', '--mapping', dest='mapping')
    parser.add_argument('-o', '--output', dest='output')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    orig_mtx = pd.read_table(args.input_matrix, index_col=0)
    mapping = pd.read_table(args.mapping, index_col=0, names=['ensg', 'symbol', 'info'], header=0)
    mapping = mapping.drop_duplicates()
    merged = pd.merge(orig_mtx, mapping, left_index=True, right_index=True)
    merged = merged.set_index('symbol')
    merged = merged.drop('info', axis=1)
    merged.to_csv(args.output, sep='\t')
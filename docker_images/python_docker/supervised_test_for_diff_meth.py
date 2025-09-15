import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, false_discovery_control


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--methylation_file', type=Path, required=True)
    parser.add_argument('-a', '--annotation_file', type=Path, required=True)
    # here this is the percentile, which depends on the `style` arg.
    # For instance, if the value here is 20 and `style` is hyper, then
    # we will get the 100-20=80th percentile.
    parser.add_argument('-s', '--style', choices=['hypo', 'hyper'], required=True)

    # what fraction of probes need to be non-NaN to keep the row
    parser.add_argument('-t', '--na_threshold', default=0.5)
    return parser.parse_args()


def run_tests(low_exp_df, high_exp_df, style):

    # save the probe IDs before we just grab the numpy matrices
    probe_ids = low_exp_df.index

    low_exp_mtx = low_exp_df.values
    high_exp_mtx = high_exp_df.values

    # we are testing whether the high group is either lower (for testing hypomethylation)
    # or greater (for testing hypermethylation)
    alternative = 'less' if style=='hypo' else 'greater'
    result = mannwhitneyu(high_exp_mtx, low_exp_mtx, alternative=alternative, axis=1, nan_policy='omit')
    pvals = result.pvalue
    fdr = false_discovery_control(pvals, method='bh')

    m1 = np.nanmean(low_exp_mtx, axis=1)
    m2 = np.nanmean(high_exp_mtx, axis=1)

    if style == 'hypo':
        effect_size = m1 - m2
    else:
        effect_size = m2 - m1

    result_df = pd.DataFrame({
        'low_group': m1,
        'high_group': m2,
        'effect_size': effect_size,
        'pval': pvals,
        'padj': fdr
    }, index=probe_ids)
    return result_df.sort_values(['padj', 'effect_size'], ascending=[True, False])


if __name__ == '__main__':

    args = parse_args()
    meth_df = pd.read_table(args.methylation_file, sep='\t', index_col=0)
    ann_df = pd.read_table(args.annotation_file, sep='\t', index_col=0)
    ann_df.index = ['-'.join(x.split('-')[:4]) for x in ann_df.index]

    # separate into samples that express high and low o the target gene
    low_exp_df = meth_df[ann_df.loc[ann_df.expression_state == 'low'].index]
    high_exp_df = meth_df[ann_df.loc[ann_df.expression_state == 'high'].index]

    # filter to remove rows with a lot of NAs, based on the `na_threshold` parameter
    N1 = low_exp_df.shape[1]
    N2 = high_exp_df.shape[1]
    low_exp_df = low_exp_df.dropna(axis=0, thresh=int(N1 * args.na_threshold))
    high_exp_df = high_exp_df.dropna(axis=0, thresh=int(N2 * args.na_threshold))

    intersection_set = (low_exp_df.index).intersection(high_exp_df.index)
    low_exp_df = low_exp_df.loc[intersection_set]
    high_exp_df = high_exp_df.loc[intersection_set]

    result = run_tests(low_exp_df, high_exp_df, args.style)

    # The names are structured like: TCGA-COAD.betas.elmer_enhancer_probes.reformatted_sample_id.tsv
    # and extracting them out lets us keep track of how the tests were run
    name_components = args.methylation_file.name.split('.')
    tcga_id = name_components[0]
    probe_class = name_components[2]
    result.to_csv(f'{tcga_id}.{probe_class}.supervised.{args.style}_na{args.na_threshold}.tsv', sep='\t')
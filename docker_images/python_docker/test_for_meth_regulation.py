import argparse
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, false_discovery_control


def parse_args():
    parser = argparse.ArgumentParser()
    # results of the differential methylation testing
    parser.add_argument('-m', '--methylation_results', type=Path, required=True)

    # the matrix of beta values
    parser.add_argument('-b', '--betas', type=Path, required=True)

    # normalized gene expression
    parser.add_argument('-n', '--norm_counts_file', type=Path, required=True)

    # a file giving the probe-to-closest-gene mappings
    parser.add_argument('-t', '--gene_targets', type=Path, required=True)

    # which samples were deemed high and low for our target gene
    parser.add_argument('-a', '--annotations', type=Path, required=True)

    parser.add_argument('-p', '--padj_threshold', type=float, default=0.01)
    parser.add_argument('-s', '--effect_size', type=float, default=0.2)

    return parser.parse_args()



if __name__ == '__main__':

    args = parse_args()
    diff_methylation_testing_df = pd.read_table(args.methylation_results, sep='\t', index_col=0)
    ann_df = pd.read_table(args.annotations, sep='\t', index_col=0)
    norm_counts = pd.read_table(args.norm_counts_file, sep='\t', index_col=0)
    target_mapping = pd.read_table(args.gene_targets, 
                                   header=None, 
                                   usecols=[3, 9],
                                   names=['probe_id','gene_id'], 
                                   sep='\t').drop_duplicates()
    betas = pd.read_table(args.betas, index_col=0, sep='\t')

    # extract the parameters for the direction of methylation. files are named like:
    # TCGA-COAD.elmer_enhancer_probes.supervised.hypo_20_na0.5.tsv
    name_contents = args.methylation_results.name.split('.')
    tcga_type = name_contents[0]
    probe_class = name_contents[1]
    supervised_status = name_contents[2]
    direction_str = name_contents[3]
    hypo_or_hyper_choice = direction_str.split('_')[0]

    # filter to keep only the probes that were differentially methylated:
    f1 = diff_methylation_testing_df.padj < args.padj_threshold
    f2 = diff_methylation_testing_df.effect_size >= args.effect_size
    diff_methylation_testing_df = diff_methylation_testing_df.loc[f1 & f2]

    if diff_methylation_testing_df.shape[0] == 0:
        sys.stdout.write('No significant probes at the thresholds chosen')
        open(f'{tcga_type}.{probe_class}.{supervised_status}.{hypo_or_hyper_choice}.probe_to_gene_results.tsv', 'w').close()
        sys.exit(0)

    # filter the beta matrix to keep only those probes of interest
    betas = betas.loc[diff_methylation_testing_df.index]

    # filter the target_mapping so we keep only the probes which are differentially methylated:
    target_mapping = target_mapping.loc[target_mapping['probe_id'].isin(diff_methylation_testing_df.index)]

    # filter the counts to only get the samples corresponding to low- and high-expression
    # of the driver gene
    nc_low = norm_counts[ann_df.index[ann_df['expression_state'] == 'low']]
    nc_high = norm_counts[ann_df.index[ann_df['expression_state'] == 'high']]

    # if the `hypo_or_hyper_choice` variable is 'hypo' then the test for differential methylation was testing
    # if the probes in the high-expressing cohort were hypomethylated. In that case,
    # we are looking to test if the nc_high values are generally higher than nc_low values. 
    # Conversely, if we are considering hypermethylated probes, we want to check if the 
    # values in nc_high are characteristically LOWER than in nc_low.
    # Note that we are performing this test for ALL the genes (not just those associated
    # with differentially expressed probes) so we can easily have this pre-processed
    # for the empirical p-value calculations
    alternative = 'greater' if hypo_or_hyper_choice=='hypo' else 'less'
    res = mannwhitneyu(nc_high, nc_low, alternative=alternative, axis=1)
    all_test_pvals = pd.DataFrame({'mw_u': res.pvalue}, index=nc_high.index)

    result_df = pd.DataFrame()
    for probe_id, subdf in target_mapping.groupby('probe_id'):
        ensg_ids = subdf.gene_id.unique()
        assoc_pvals = all_test_pvals.loc[ensg_ids, ['mw_u']].copy()
        assoc_pvals['probe_id'] = probe_id
        assoc_pvals = assoc_pvals.reset_index(names='ensg_id')
        result_df = pd.concat([result_df, assoc_pvals], axis=0)
    result_df['fdr'] = false_discovery_control(result_df.mw_u, method='bh')
    result_df.to_csv(f'{tcga_type}.{probe_class}.{supervised_status}.{hypo_or_hyper_choice}.probe_to_gene_results.tsv', sep='\t', index=False)

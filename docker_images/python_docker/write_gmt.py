import argparse
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mapping',
                      type=Path, required=True)
    parser.add_argument('infiles', 
                        nargs='+')
    parser.add_argument('-p', '--padj_threshold', type=float, default=0.001)

    return parser.parse_args()


def write_gmt_file(gene_sets, output_fname):
    '''
    Write a GMT-format file (e.g. like msigDB)
    `gene_sets` is a list of sets.
    '''
    with open(output_fname, 'w') as fout:
        for k, gene_set in gene_sets.items():
            gs = '\t'.join(gene_set)
            fout.write(f'{k}\t-\t{gs}\n')


if __name__ == '__main__':

    args = parse_args()

    mapping = pd.read_table(args.mapping, sep='\t')


    gene_sets = {}
    tcga_types = set()
    probe_styles = set()
    for f in args.infiles:
        fname = f.split('/')[-1]
        tcga_type, probe_style, _a, direction, _b, _c = fname.split('.')
        tcga_types.add(tcga_type)
        probe_styles.add(probe_style)
        df = pd.read_table(f)
        df = df.loc[df.fdr <= args.padj_threshold]

        # the GO analysis will require gene symbols, so we map that here
        df = pd.merge(df, mapping, left_on='ensg_id', right_on='gene_id')
        
        gene_set_name = f'gs_{direction}'
        gene_sets[gene_set_name] = set(df['gene_name'])

    assert(len(tcga_types) == 1)
    assert(len(probe_styles) == 1)

    output_name = f'{list(tcga_types)[0]}.{list(probe_styles)[0]}.gmt'
    write_gmt_file(gene_sets, output_name)


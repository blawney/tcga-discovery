import argparse
from pathlib import Path
import pickle

import pandas as pd
import numpy as np
import networkx as nx


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--tcga_type', 
                        required=True)
    parser.add_argument('--interaction_matrix', 
                        required=True,
                        type=Path)
    parser.add_argument('--info_matrix',
                        required=True,
                        type=Path)
    parser.add_argument('--dge_matrix',
                        required=True,
                        type=Path)
    parser.add_argument('--score_threshold',
                        default=900,
                        type=int)
    parser.add_argument('--padj_threshold',
                        default=0.05,
                        type=float)
    parser.add_argument('--min_community_size',
                        default=15,
                        type=int)
    return parser.parse_args()


def has_both_proteins(row, pset):
    has_a = row.protein1 in pset
    has_b = row.protein2 in pset
    return has_a & has_b


def has_at_least_one_protein(row, pset):
    has_a = row.protein1 in pset
    has_b = row.protein2 in pset
    return has_a or has_b


def write_gmt_file(gene_sets, output_fname):
    '''
    Write a GMT-format file (e.g. like msigDB)
    `gene_sets` is a list of sets.
    '''
    with open(output_fname, 'w') as fout:
        for i, gene_set in enumerate(gene_sets):
            gs = '\t'.join(gene_set)
            fout.write(f'gs{i}\t-\t{gs}\n')

def create_graph(interaction_df, tcga_type, graph_label):
    # create the interaction graph
    G = nx.Graph()
    G.add_edges_from(interaction_df[['preferred_name1','preferred_name2']].to_records(index=False))

    # find communities and write to a GMT-format file:
    communities = nx.community.louvain_communities(G)
    communities = [x for x in communities if len(x) >= args.min_community_size]
    write_gmt_file(communities, f'{tcga_type}.louvain_communities.{graph_label}.gmt')

    # collect some basic node metrics and save
    btn_c = pd.Series(nx.betweenness_centrality(G), name='betweenness')
    degree_c = pd.Series(nx.degree_centrality(G), name='degree')
    evector_c = pd.Series(nx.eigenvector_centrality(G, max_iter=1000), name='evector_centrality')
    stats_df = pd.concat([btn_c, degree_c, evector_c], axis=1)
    stats_df.to_csv(f'{tcga_type}.network_stats.{graph_label}.tsv', sep='\t')
    
    # save the graph as a pickle
    with open(f'{tcga_type}.graph.{graph_label}.pkl', "wb") as f:
        pickle.dump(G, f)


if __name__ == '__main__':

    args = parse_args()

    # this has the two ENSP proteins and their interaction score
    full_interaction_df = pd.read_table(args.interaction_matrix, sep=' ')
    full_interaction_df = full_interaction_df.loc[full_interaction_df.combined_score > args.score_threshold]

    # this has the mapping from ENSP (#string_protein_id) to 
    # common gene symbol (preferred_name)
    mapping_df = pd.read_table(args.info_matrix, usecols=[0,1], index_col=0)

    # differential expression results
    dge_df = pd.read_table(args.dge_matrix, index_col=0)
    dge_df = dge_df.loc[dge_df.padj < args.padj_threshold]

    # by merging the dge matrix and mapping, we keep only gene/protein symbols
    # corresponding to differentially expressed genes
    mdf = pd.merge(dge_df, mapping_df, left_index=True, right_on='preferred_name')
    diff_protein_set = set(mdf.index.unique())

    interaction_filter = full_interaction_df.apply(has_at_least_one_protein, axis=1, args=(diff_protein_set,))
    interaction_df = full_interaction_df.loc[interaction_filter]

    # map both proteins to the gene symbols:
    interaction_df = pd.merge(interaction_df, mapping_df, left_on='protein1', right_index=True)
    interaction_df = pd.merge(interaction_df, mapping_df, left_on='protein2', right_index=True, suffixes=['1','2'])

    create_graph(interaction_df, args.tcga_type, 'full')

    interaction_filter = full_interaction_df.apply(has_both_proteins, axis=1, args=(diff_protein_set,))
    interaction_df = full_interaction_df.loc[interaction_filter]

    # map both proteins to the gene symbols:
    interaction_df = pd.merge(interaction_df, mapping_df, left_on='protein1', right_index=True)
    interaction_df = pd.merge(interaction_df, mapping_df, left_on='protein2', right_index=True, suffixes=['1','2'])

    create_graph(interaction_df, args.tcga_type, 'dge_only')
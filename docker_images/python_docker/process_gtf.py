import argparse

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_gtf', required=True)
    parser.add_argument('-b', '--bed_output', required=True)
    parser.add_argument('-m', '--map_output', required=True)
    return parser.parse_args()


def extract_field_from_gtf(info, fields_to_extract):
    contents = [x.strip() for x in info.split(';')]
    contents = [x for x in contents if len(x) > 0]
    d = {x:y[1:-1] for x,y in (z.split(' ') for z in contents)}
    return pd.Series({k:d[k] for k in fields_to_extract})


if __name__ == '__main__':
    args = parse_args()

    gtf_df = pd.read_table(args.input_gtf, 
                           sep='\t', 
                           comment='#', 
                           low_memory=False,
                           names=['chrom','source','feature','start','end','_0','strand', '_1', 'info'])
    
    # filter to keep only gene features:
    gtf_df = gtf_df.loc[gtf_df.feature == 'gene']

    gene_info = gtf_df['info'].apply(lambda x: extract_field_from_gtf(x, ['gene_id', 'gene_name']))
    gene_info['gene_id'] = gene_info['gene_id'].apply(lambda x: x.split('.')[0])

    gtf_df['ensg_id'] = gene_info['gene_id']

    # now save the BED file of genic regions:
    gtf_df[['chrom','start','end', 'ensg_id', '_0', 'strand']].to_csv(args.bed_output, 
                                                                      header=False, 
                                                                      sep='\t', 
                                                                      index=False)

    # save the ENSG to symbol mapping
    gene_info[['gene_id', 'gene_name']].to_csv(args.map_output, sep='\t', index=False)
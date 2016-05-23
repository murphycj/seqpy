"""
This script combines the gene and isoform fpkm estimates from cufflinks
from different samples into one file

"""

import sys
import pandas
import numpy as np
import os
import argparse

def main(args):

    data = zip(args.files,args.samples)

    results = {}
    genes = []
    samples = []

    for d in data:

        fpkm_data = pandas.read_table(d[0], sep='\t', index_col=0)


        if args.remove_duplicates:
            genes += fpkm_data.index.tolist()
            temp2 = pandas.DataFrame(fpkm_data['FPKM'])
            temp2['index'] = temp2.index
            temp2.drop_duplicates(subset='index', take_last=True, inplace=True)
            del temp2['index']
            results[d[1]] = temp2
        else:
            pass

        samples.append(d[1])

    genes = list(set(genes))
    for s, fpkm_data in results.items():
        results[s] = fpkm_data.ix[genes]

    all_fpkm_data = pandas.DataFrame(index = list(set(genes)), columns = samples)
    all_fpkm_data.values.fill(0.0)
    for s, fpkm_data in results.items():
        all_fpkm_data[s] = fpkm_data
    all_fpkm_data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument(
    '--files',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of cufflinks output files.'
)
parser.add_argument(
    '--samples',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of sample names'
)
parser.add_argument(
    '--remove_duplicates',
    action='store_true',
    help='Remove duplicate gene names (default True), ' + \
         'otherwise keep duplicate gene names by indexing them (e.g.' + \
         'TP53-1,TP53-2,...)',
    default=True,
    required=False
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

main(args=args)

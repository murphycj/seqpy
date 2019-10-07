"""
This script combines the gene and isoform fpkm estimates from cufflinks
from different samples into one file

"""

import sys
import os
import argparse

import pandas as pd
import numpy as np

def main(args):

    data = zip(args.files,args.samples)
    results = {}

    for d in data:

        fpkm_data = pd.read_csv(d[0], sep='\t')
        temp2 = pd.DataFrame(fpkm_data[['tracking_id','FPKM']])

        if args.duplicates=='random':
            temp2.drop_duplicates(subset='tracking_id', take_last=True, inplace=True)
            temp2.index = temp2['tracking_id']
            results[d[1]] = temp2['FPKM']
        elif args.duplicates=='highest':
            temp2 = temp2.groupby('tracking_id',group_keys=False).apply(lambda x: x.ix[x.FPKM.idxmax()])
            temp2.index = temp2['tracking_id']
            results[d[1]] = temp2['FPKM']
        else:
            print 'Nothing done!'
            sys.exit()

    for s, fpkm_data in results.items():
        results[s] = fpkm_data.to_dict()

    all_fpkm_data = pd.DataFrame(results)
    all_fpkm_data.fillna(0.0)
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
    '--duplicates',
    type=str,
    help='How to handle duplicates: random (randomly select one), index (e.g. TP53-1,TP53-2,...), or highest (default, take one with highest FPKM).',
    default='highest',
    required=False
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

if args.duplicates not in ['random','highest','index']:
    print 'Select one of: random, highest, or index for --duplicates.'
    sys.exit()

main(args=args)

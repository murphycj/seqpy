"""
This will aggregate the output from HTSeq-count on multiple samples
into one file

"""

import sys
import pandas as pd
import numpy as np
import os
import argparse

def main(args):
    data = pd.DataFrame()

    for fpath, name in zip(args.files,args.samples):

        sample_data = pd.read_csv(fpath, sep='\t',header=None)
        sample_data = sample_data[~sample_data[0].str.contains('__')]
        sample_data.index = sample_data[0]
        sample_data = sample_data[[1]]
        sample_data.columns = [name]

        data = pd.concat([data, sample_data], axis=1)

    data = data.fillna(0)
    data.index.name = None
    
    data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Aggregates output from HTSeq-count')
parser.add_argument(
    '--files',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of output files.'
)
parser.add_argument(
    '--samples',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of sample names'
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

main(args=args)

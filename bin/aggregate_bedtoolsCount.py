"""
This will aggregate the output from HTSeq-count on multiple samples
into one file

"""

import sys
import pandas
import numpy as np
import os
import argparse

def main(args):

    data = zip(args.files,args.samples)

    sample_data = pandas.read_table(
        data[0][0],
        sep='\t',
        index_col=None,
        header=None,
        low_memory=False
    )

    all_data = pandas.DataFrame(
        index=range(0,sample_data.shape[0]),
        columns=['chr','start','end'] + args.samples
    )
    all_data['chr'] = sample_data[0]
    all_data['start'] = sample_data[1]
    all_data['end'] = sample_data[2]

    for d in data:

        sample_data = pandas.read_table(
            d[0],
            sep='\t',
            index_col=None,
            header=None,
            low_memory=False
        )

        all_data[d[1]] =  sample_data[3]
    all_data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Aggregates output from Bedtools count')
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

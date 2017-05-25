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

    results = {}
    genes = set()
    samples = []

    for d in zip(args.files,args.samples):

        sample_data = pandas.read_table(d[0], sep='\t',header=None)
        sample_data = sample_data[~sample_data[0].str.contains('__')]
        genes.update(sample_data[0].tolist())
        sample_data.index = sample_data[0]
        sample_data = sample_data[1]
        results[d[1]] = sample_data
        samples.append(d[1])

    genes = list(genes)

    all_data = pandas.DataFrame(index = list(set(genes)), columns = samples)
    all_data.values.fill(0.0)
    for s, d in results.items():
        all_data[s] = d.ix[genes]
    all_data.to_csv(args.out)


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

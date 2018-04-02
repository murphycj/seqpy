
import sys
import pandas
import numpy as np
import os
import argparse

def main(args):

    data = zip(args.dirs,args.samples)

    results = {}
    genes = []
    samples = []

    for d in data:

        tpm_data = pandas.read_table(d[0] + '/abundance.tsv', sep='\t')
        tpm_data.index = tpm_data['target_id']
        genes += tpm_data['target_id'].tolist()
        temp2 = pandas.DataFrame(tpm_data[['target_id','tpm']])
        results[d[1]] = temp2['tpm']
        samples.append(d[1])

    genes = list(set(genes))

    for s, tpm_data in results.items():
        results[s] = tpm_data.ix[genes]

    all_tpm_data = pandas.DataFrame(index = list(set(genes)), columns = samples)
    all_tpm_data.values.fill(0.0)
    for s, tpm_data in results.items():
        all_tpm_data[s] = tpm_data
    all_tpm_data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Aggregates the tpm values from kallisto output')
parser.add_argument(
    '--dirs',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of cufflinks output dirs.'
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

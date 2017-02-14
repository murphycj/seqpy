"""
This script combines the gene and isoform fpkm estimates from cufflinks
from different samples into one file

"""

import glob
import sys
import pandas
import os
import argparse

suffixes = [
    'sample_cumulative_coverage_counts',
    'sample_cumulative_coverage_proportions',
    'sample_interval_statistics',
    'sample_interval_summary',
    'sample_statistics',
    'sample_summary'
    ]

def aggregate_summary(args):

    inputs = zip(args.sample_dir,args.samples)

    results = {}
    for suffix in suffixes:
        results[suffix] = {}

        for dir_name, sample in inputs:

            fname = glob.glob(os.path.join(dir_name,"*" + suffix))
            assert len(fname)==1, 'Too many or too few files found: ' + str(fname)

            data = pandas.read_table(fname[0], sep='\t', index_col=0)

            if suffix=='sample_statistics':

                results[suffix][sample] = data.ix['sample_test']
            elif suffix=='sample_interval_summary':
                results[suffix][sample] = data
            elif suffix=='sample_interval_statistics':
                results[suffix][sample] = data.ix['At_least_1_samples']
            elif suffix=='sample_cumulative_coverage_counts':
                results[suffix][sample] = data.ix['NSamples_1']
            else:
                try:
                    results[suffix][sample] = data.ix[data.index[0]]
                except:
                    import pdb; pdb.set_trace()


    all_data = pandas.DataFrame(results['sample_statistics'])
    all_data.to_csv(args.prefix + '.sample_statistics.csv')

    all_data = pandas.DataFrame(results['sample_summary'])
    all_data.to_csv(args.prefix + '.sample_summary.csv')

    all_data = pandas.DataFrame(results['sample_interval_statistics'])
    all_data.to_csv(args.prefix + '.sample_interval_statistics.csv')

    all_data = pandas.DataFrame(results['sample_cumulative_coverage_proportions'])
    all_data.to_csv(args.prefix + '.sample_cumulative_coverage_proportions.csv')

    all_data = pandas.DataFrame(results['sample_cumulative_coverage_counts'])
    all_data.to_csv(args.prefix + '.sample_cumulative_coverage_counts.csv')

    for c in ['_total_cvg','_mean_cvg','_granular_Q1','_granular_median','_granular_Q3']:
        r = {}
        for sample, data in results['sample_interval_summary'].items():

            r[sample] = data[sample + c]

        all_data = pandas.DataFrame(r)
        all_data.to_csv(args.prefix + '.sample_interval_summary.' + c + '.csv')

parser = argparse.ArgumentParser(description='Aggregate summary files')
parser.add_argument(
    '--sample_dir',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of directories containing GATK coverage output'
)
parser.add_argument(
    '--samples',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of sample names'
)
parser.add_argument(
    '--prefix',
    type=str,
    required=True,
    help='Output prefix'
)

args = parser.parse_args()

aggregate_summary(args=args)

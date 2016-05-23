"""
Combine the "Log.final.out" files from STAR mapping into one CSV file

"""

import argparse
import pandas
from crimson import flagstat
import glob
import os

VARIABLES = [
    'Number of input reads','Average input read length','Uniquely mapped reads number',
    'Uniquely mapped reads %','Average mapped length','Number of splices: Total',
    'Number of splices: Annotated (sjdb)','Number of splices: GT/AG',
    'Number of splices: GC/AG','Number of splices: AT/AC',
    'Number of splices: Non-canonical','Mismatch rate per base, %',
    'Deletion rate per base','Deletion average length','Insertion rate per base',
    'Insertion average length','Number of reads mapped to multiple loci',
    '% of reads mapped to multiple loci','Number of reads mapped to too many loci',
    '% of reads mapped to too many loci','% of reads unmapped: too many mismatches',
    '% of reads unmapped: too short','% of reads unmapped: other'
]

def parse_file(filename):
    """
    Parse the individual file getting the data
    """

    data = dict(zip(VARIABLES,[0]*len(VARIABLES)))

    fin = open(filename,'r')
    for line in fin.readlines():
        line = ' '.join(line.rstrip().split())
        for i in VARIABLES:
            if line.find(i)!=-1:
                try:
                    data[i] = line.split(' | ')[1]
                except:
                    import pdb; pdb.set_trace()

    return data

def main(args):
    data = zip(args.files,args.samples)
    samples = []
    results = dict()

    for d in data:
        results[d[1]] = parse_file(d[0])
        samples.append(d[1])

    all_data = pandas.DataFrame(
        index = VARIABLES,
        columns = samples
    )

    results = pandas.DataFrame(results)
    results = results.ix[VARIABLES]
    results.to_csv(args.out)

parser = argparse.ArgumentParser(description='Combine STAR Log.final.out files')
parser.add_argument(
    '--files',
    nargs='+',
    type=str,
    required=True,
    help='Space-delimited list of files'
)
parser.add_argument(
    '--samples',
    nargs='+',
    type=str,
    required=True,
    help='Space-delimited list of files'
)
parser.add_argument(
    '--out',
    type=str,
    help='Output file name',
    required=True
)
args = parser.parse_args()

main(args=args)

"""
Aggregates flagstat output files from samtools into a single csv file

"""


import argparse
import pandas
from crimson import flagstat
import glob
import os

def main(args):
    data = {}

    for f in args.flagstat:
        sample = os.path.split(f)[1]
        data[sample] = flagstat.parse(f)['pass_qc']
    results = pandas.DataFrame(data).T
    results.to_csv(args.out)

parser = argparse.ArgumentParser(description='Combine flagstat files')
parser.add_argument('--flagstat',nargs='+',help='Space-delimited list of flagstat files')
parser.add_argument('--out',type=str,help='Output file name',required=True)
args = parser.parse_args()

main(args=args)

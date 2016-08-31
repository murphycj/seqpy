import sys
import pandas
import os
import argparse

def main(args):

    counts = []

    data = zip(args.counts,args.samples)

    samples = []

    for c,s in data:
        temp = pandas.read_csv(c,index_col=0)
        samples.append(s)
        counts.append(temp)

    data = pandas.concat(counts,axis=1)
    data.columns = samples

    data.to_csv(args.out)

parser = argparse.ArgumentParser(description='')
parser.add_argument('--counts',type=str,required=True,help='Space-delimited list of count files',nargs='+')
parser.add_argument('--samples',type=str,required=True,help='Sample names (same order as counts)',nargs='+')
parser.add_argument('--out',type=str,help='out file',required=False)
args = parser.parse_args()

main(args=args)

import sys
import pandas
import numpy as np
import glob
import os
import argparse

def get_htseqcount(args):

    results = {}
    genes = []
    samples = []

    dirs = glob.glob(args.dir + '/' + args.dirPattern + '*')
    for f in dirs:
        if not os.path.isdir(f):
            continue
        sample = os.path.split(f)[1]

        if sample in args.skip:
            continue

        files = glob.glob(f + '/*' + args.countPattern)

        if len(files)==0:
            print "No file for " + sample
            continue

        if len(files) > 1:
            print "More than one file found: " + sample
            sys.exit(1)
        try:
            d = pandas.read_table(files[0], sep='\t', index_col=0,header=None)
        except:
            print sample
            continue


        d = d.ix[d.index.map(lambda x: x.find('__')==-1)]
        genes += d.index.tolist()
        results[sample] = d
        samples.append(sample)

    data = pandas.DataFrame(index = list(set(genes)), columns = samples)
    data.values.fill(0.0)
    for s, d in results.items():
        data[s] = d
    data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Aggregates count data')
parser.add_argument('--dir',type=str,help='Meta-directory where each sub-directory contains fpkm for each sample',required=True)
parser.add_argument('--dirPattern',type=str,help='Directory pattern',required=True)
parser.add_argument('--countPattern',type=str,help='Pattern to use when searching for the count files',required=True)
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--skip',nargs='*',help='Specific samples to skip',required=False,default=[])
args = parser.parse_args()


get_htseqcount(args=args)

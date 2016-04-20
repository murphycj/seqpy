import sys
import pandas
import numpy as np
import glob
import os
import argparse

def parse_metrics_file(filename):

    results = {}

    parse_next = False

    fin = open(filename,'r')
    for line in fin.readlines():
        if 'PF_BASES' in line:
            temp = line.split('\t')
            parse_next=True
        elif parse_next:
            parse_next=False


    return results


def main(args):

    results = {}
    genes = []
    samples = []

    dirs = glob.glob(args.dir + '/*')

    for f in dirs:
        if not os.path.isdir(f):
            continue

        sample = os.path.split(f)[1]

        if sample in args.skip:
            continue

        files = ''
        if not args.walk:

            files = glob.glob(f + '/*RNAseq_metrics.txt')

            if len(files)==0:
                print "No file for " + sample
                continue

            if len(files) > 1:
                print "More than one file found: " + sample
                sys.exit(1)

            files = files[0]
        else:
            pattern = 'RNAseq_metrics.txt'

            nf = 0
            for root, subdir, filename in os.walk(f):
                for f1 in filename:
                    if pattern in f1:
                        nf += 1
                        files = os.path.join(root, f1)
            if nf > 1:
                print "More than one file found: " + sample
                sys.exit(1)

        if files == '':
            print 'No file found'
            continue

        fin = open(files,'r').read().split('\n')
        get_next = False
        columns = []
        values = []
        for i in range(0,len(fin)):
            if 'PF_BASES' in fin[i]:
                columns = fin[i].rstrip().split('\t')
                values = fin[i+1].rstrip().split('\t')
                break
        results[sample] = dict(zip(columns,values))

    data = pandas.DataFrame(results).T
    data.to_csv(args.out)



parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument('--dir',type=str,help='Meta-directory where each sub-directory contains fpkm for each sample',required=True)
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--skip',nargs='*',help='Specific samples to skip',required=False,default=[])
parser.add_argument('--walk',action='store_true',help='Recursively walk through the sample dirs to find the file (defaule false).',required=False)
args = parser.parse_args()

main(args=args)

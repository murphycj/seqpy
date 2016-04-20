from crimson import fastqc
import sys
import pandas
import numpy as np
import glob
import os
import argparse

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

            files = glob.glob(f + '/*stdin_fastqc.zip')

            if len(files)==0:
                print "No file for " + sample
                continue

            if len(files) > 1:
                print "More than one file found: " + sample
                sys.exit(1)

            files = files[0]
        else:
            pattern = 'stdin_fastqc.zip'

            nf = 0
            for root, subdir, file in os.walk(f):
                if pattern in file:
                    nf += 1
                    files = os.path.join(root, pattern)
            if nf > 1:
                print "More than one file found: " + sample
                sys.exit(1)
        if files == '':
            continue
        os.system('unzip -o -q ' + os.path.abspath(files) + ' -d tmp-fastqc')
        f = fastqc.parse('./tmp-fastqc/stdin_fastqc/fastqc_data.txt')

        results[sample] = {
            'Per base sequence quality':f['Per base sequence quality']['status'],
            'Sequence Duplication Levels':f['Sequence Duplication Levels']['status'],
            'Per sequence GC content':f['Per sequence GC content']['status'],
            'Sequence Length Distribution':f['Sequence Length Distribution']['status'],
            'Kmer Content':f['Kmer Content']['status'],
            'Overrepresented sequences':f['Overrepresented sequences']['status'],
            'Per base N content':f['Per base N content']['status'],
            'Per sequence quality scores':f['Per sequence quality scores']['status']
        }
    data = pandas.DataFrame(results)
    data.to_csv(args.out)

    os.system('rm -rf tmp-fastqc')


parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument('--dir',type=str,help='Meta-directory where each sub-directory contains fpkm for each sample',required=True)
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--skip',nargs='*',help='Specific samples to skip',required=False,default=[])
parser.add_argument('--walk',action='store_true',help='Recursively walk through the sample dirs to find the file (defaule false).',required=False)
args = parser.parse_args()

main(args=args)

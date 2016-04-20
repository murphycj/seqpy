import sys
import pandas
import numpy as np
import glob
import os
import argparse

def get_fpkm(args):

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
            if args.genes:
                files = glob.glob(f + '/genes.fpkm_tracking')
            else:
                files = glob.glob(f + '/isoforms.fpkm_tracking')

            if len(files)==0:
                print "No file for " + sample
                continue

            if len(files) > 1:
                print "More than one file found: " + sample
                sys.exit(1)

            files = files[0]
        else:
            pattern = ''
            if args.genes:
                pattern = 'genes.fpkm_tracking'
            else:
                pattern = 'isoforms.fpkm_tracking'

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

        try:
            d = pandas.read_table(files, sep='\t', index_col=0)
        except:
            print 'Could not read table'
            import pdb; pdb.set_trace()

        genes += d.index.tolist()
        temp2 = pandas.DataFrame(d['FPKM'])
        temp2['index'] = temp2.index
        temp2.drop_duplicates(subset='index', take_last=True, inplace=True)
        del temp2['index']
        results[sample] = temp2
        samples.append(sample)

    genes = list(set(genes))
    for s, d in results.items():
        results[s] = d.ix[genes]

    data = pandas.DataFrame(index = list(set(genes)), columns = samples)
    data.values.fill(0.0)
    for s, d in results.items():
        data[s] = d
    data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument('--dir',type=str,help='Meta-directory where each sub-directory contains fpkm for each sample',required=True)
parser.add_argument('--genes',action='store_true',help='Collate FPKM values for genes',required=False)
parser.add_argument('--isoforms',action='store_true',help='Collate FPKM values for isoforms',required=False)
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--skip',nargs='*',help='Specific samples to skip',required=False,default=[])
parser.add_argument('--walk',action='store_true',help='Recursively walk through the sample dirs to find the file (defaule false).',required=False)
args = parser.parse_args()

if not args.genes and not args.isoforms:
    print "Specify either --genes or --isoforms"
    sys.exit()

if args.genes and args.isoforms:
    print "Specify only --genes or --isoforms"
    sys.exit()

get_fpkm(args=args)

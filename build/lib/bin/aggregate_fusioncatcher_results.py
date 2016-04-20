import sys
import pandas
import numpy as np
import glob
import os
import argparse

def get_fusioncatcher(args):

    results = {}
    genes = []
    samples = []
    fusions = {}
    all_fusions = []
    all_data = pandas.DataFrame()

    dirs = glob.glob(args.dir + '/*')

    for f in dirs:
        if not os.path.isdir(f):
            continue

        sample = os.path.split(f)[1].split('-')[0]

        if sample in args.skip:
            continue

        files = ''
        if not args.walk:
            files = glob.glob(f + '/final-list_candidate-fusion-genes.txt')
            if len(files)==0:
                print "No file for " + sample
                continue
            if len(files) > 1:
                print "More than one file found: " + sample
                sys.exit(1)
            files=files[0]
        else:
            nf = 0
            for root, subdir, file in os.walk(f):
                if 'final-list_candidate-fusion-genes.txt' in file:
                    nf += 1
                    files = os.path.join(root, 'final-list_candidate-fusion-genes.txt')
            if nf > 1:
                print "More than one file found: " + sample
                sys.exit(1)
            elif nf==0:
                continue

        fusions[sample] = []

        data = pandas.read_table(files, sep='\t')

        col_data = []
        for i in data.index:
            gene1 = data.ix[i]['Gene_1_symbol(5end_fusion_partner)']
            gene2 = data.ix[i]['Gene_2_symbol(3end_fusion_partner)']
            fusion = gene1 + '-' + gene2
            fusions[sample].append(fusion)

            if fusion not in all_fusions:
                all_fusions.append(fusion)

        data['sample'] = sample
        all_data = all_data.append(data)

    #reorder column names

    temp = all_data.columns.tolist()
    del temp[-1]
    temp = ['sample'] + temp
    all_data = all_data[temp]
    all_data.to_csv(args.out,index=False)


parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument('--dir',type=str,help='Meta-directory where each sub-directory contains fpkm for each sample',required=True)
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--skip',nargs='*',help='Specific samples to skip',required=False,default=[])
parser.add_argument('--walk',action='store_true',help='Recursively walk through the sample dirs to find the file (defaule false).',required=False)
args = parser.parse_args()

get_fusioncatcher(args=args)

"""
Aggregates the output files from fusioncatcher on different samples
into one file

"""

import sys
import pandas
import numpy as np
import glob
import os
import argparse

def main(args):

    data = zip(args.files,args.samples)

    fusions = {}
    all_fusions = []
    all_fusioncatcher_data = pandas.DataFrame()

    for d in data:

        fusions[data[0]] = []

        fusioncatcher_data = pandas.read_table(data[1], sep='\t')

        col_data = []
        for i in data.index:
            gene1 = fusioncatcher_data.ix[i]['Gene_1_symbol(5end_fusion_partner)']
            gene2 = fusioncatcher_data.ix[i]['Gene_2_symbol(3end_fusion_partner)']
            fusion = gene1 + '-' + gene2
            fusions[sample].append(fusion)

            if fusion not in all_fusions:
                all_fusions.append(fusion)

        fusioncatcher_data['sample'] = sample
        all_fusioncatcher_data = all_data.append(fusioncatcher_data)

    #reorder column names

    temp = all_fusioncatcher_data.columns.tolist()
    del temp[-1]
    temp = ['sample'] + temp
    all_fusioncatcher_data = all_fusioncatcher_data[temp]
    all_fusioncatcher_data.to_csv(args.out,index=False)


parser = argparse.ArgumentParser(description='Aggregates fusioncatcher output from different samples into one file')
parser.add_argument('--files',type=str,required=True,nargs='+',help='Space-delimited list of fusioncatcher output files.')
parser.add_argument('--samples',type=str,required=True,nargs='+',help='Space-delimited list of sample names')
parser.add_argument('--out',type=str,help='Output file name',required=True)
args = parser.parse_args()

main(args=args)

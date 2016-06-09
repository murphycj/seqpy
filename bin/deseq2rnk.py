"""
Converts the CSV output file from DESeq2 or limma voom into an
RNK file for use in GSEA preranked tool
"""

import numpy as np
import argparse
import os
import pandas
import mouse2human

def create_rnk(args,data):

    result_data = pandas.DataFrame(columns=['name','rank'],index=data.index)

    if args.limmavoom:
        if args.ranking == 'statistic':
            result_data['name'] = data.index
            result_data['rank'] = data['t'].tolist()
            result_data = result_data[
                (~pandas.isnull(result_data['name'])) &
                (~pandas.isnull(result_data['rank'])) &
                (~np.isinf(result_data['rank']))
                ]
        elif args.ranking == 'logpvalFC':
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['P.Value'])*data['logFC']).tolist()
            result_data = result_data[
                (~pandas.isnull(result_data['name'])) &
                (~pandas.isnull(result_data['rank'])) &
                (~np.isinf(result_data['rank']))
                ]
    else:
        if args.ranking == 'statistic':
            result_data['name'] = data.index
            result_data['rank'] = data['stat'].tolist()
            result_data = result_data[
                (~pandas.isnull(result_data['name'])) &
                (~pandas.isnull(result_data['rank'])) &
                (~np.isinf(result_data['rank']))
                ]
        elif args.ranking == 'logpvalFC':
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['pvalue'])*data['log2FoldChange']).tolist()
            result_data = result_data[
                (~pandas.isnull(result_data['name'])) &
                (~pandas.isnull(result_data['rank'])) &
                (~np.isinf(result_data['rank']))
                ]

    result_data.to_csv(args.out,index=False,header=False, sep='\t')

def main(args):

    data = pandas.read_csv(args.infile, header=0, index_col=0)
    data = data.dropna()

    create_rnk(args=args,data=data)

parser = argparse.ArgumentParser(description='Creates RNK file from DESeq2 or LimmaVoom output.')
parser.add_argument('--infile',type=str,help='csv file containing deseq2 results',required=True)
parser.add_argument('--ranking',type=str,help='Ranking statistic to use (statistic (default), logpvalFC)',default='statistic')
parser.add_argument('--limmavoom',action='store_true',help='Include if results are from limmavoom')
parser.add_argument('--out',type=str,help='Output prefix',required=True)
args = parser.parse_args()

main(args=args)

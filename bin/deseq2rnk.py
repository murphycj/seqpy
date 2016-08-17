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
        elif args.ranking == 'logPvalSignFC':
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['P.Value'])*np.sign(data['logFC'])).tolist()
            result_data = result_data[
                (~pandas.isnull(result_data['name'])) &
                (~pandas.isnull(result_data['rank'])) &
                (~np.isinf(result_data['rank']))
                ]
        else:
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['P.Value'])).tolist()
            result_data = result_data[
                (~pandas.isnull(result_data['name'])) &
                (~pandas.isnull(result_data['rank'])) &
                (~np.isinf(result_data['rank']))
                ]

    else:
        if args.ranking == 'statistic':
            result_data['name'] = data.index
            result_data['rank'] = data['stat'].tolist()
        elif args.ranking == 'logPvalSignFC':
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['pvalue'])*np.sign(data['log2FoldChange'])).tolist()
        elif args.ranking == 'log2FC':
            result_data['name'] = data.index
            result_data['rank'] = data['log2FoldChange'].tolist()
        else:
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['pvalue'])).tolist()


    if args.reverseSign:
        result_data['rank'] = -1 * result_data['rank']

    if args.abs:
        result_data['rank'] = abs(result_data['rank'])

    result_data.to_csv(args.out,index=False,header=False, sep='\t')

def main(args):

    data = pandas.read_csv(args.infile, header=0, index_col=0)

    if args.limmavoom:
        data['P.Value']=data['P.Value'].fillna(1.0)
        data['P.Value'] = data['P.Value'].replace(0,min(data[data['P.Value']>0]['P.Value']))

        data['padj']=data['padj'].fillna(1.0)
        data['padj'] = data['padj'].replace(0,min(data[data['padj']>0]['padj']))

        data['logFC']=data['logFC'].fillna(0.0)
        data['t']=data['t'].fillna(0.0)
    else:
        data['pvalue']=data['pvalue'].fillna(1.0)
        data['pvalue'] = data['pvalue'].replace(0,min(data[data['pvalue']>0]['pvalue']))

        data['padj']=data['padj'].fillna(1.0)
        data['padj'] = data['padj'].replace(0,min(data[data['padj']>0]['padj']))

        data['log2FoldChange']=data['log2FoldChange'].fillna(0.0)
        data['stat']=data['stat'].fillna(0.0)

    create_rnk(args=args,data=data)

parser = argparse.ArgumentParser(description='Creates RNK file from DESeq2 or LimmaVoom output.')
parser.add_argument('--infile',type=str,help='csv file containing deseq2 results',required=True)
parser.add_argument('--ranking',type=str,help='Ranking statistic to use (statistic (default), logPvalSignFC, log2FC, logPval)',default='statistic',required=False)
parser.add_argument('--limmavoom',action='store_true',help='Include if results are from limmavoom (default false)',required=False,default=False)
parser.add_argument('--abs',action='store_true',help='Take absolute value of ranking statistic (default false)',required=False,default=False)
parser.add_argument('--reverseSign',action='store_true',help='Reverse the sign of the ranking metric',required=False,default=False)
parser.add_argument('--out',type=str,help='Output prefix',required=True)
args = parser.parse_args()

main(args=args)

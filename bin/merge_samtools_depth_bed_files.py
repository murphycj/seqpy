"""
This was written to combined the bed output files from running
samtools depth.

"""

import sys
import pandas
import os
import argparse

def main(args):

    results = {}
    n = 0

    for f in args.files:

        results[n] = pandas.read_table(f,sep='\t',header=None)
        n+=1


    #merge the results

    data = pandas.merge(results[0],results[1],on=[0,1])
    for i in results.keys():
        if i==0 or i==1:
            continue
        data = pandas.merge(data,results[i],on=[0,1])

    data.to_csv(args.out,index=None,header=None,sep='\t')

parser = argparse.ArgumentParser(description='Aggregates samtools depth bed files. Assumes all bed files query the same positions.')
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--files',type=str,required=True,nargs='+')
args = parser.parse_args()

main(args=args)

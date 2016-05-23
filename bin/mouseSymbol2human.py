"""
Convert mouse gene symbols to human gene symbols
"""

import numpy as np
import argparse
import os
import pandas
import mouse2human
import argparse

def main(args):

    data = pandas.read_csv(args.infile, header=0, index_col=0)

    c = mouse2human.Conversion()

    data = c.convert_mouse2human(data=data,remove_duplicates=True)

    data.to_csv(args.out)

parser = argparse.ArgumentParser(description='Converts mouse genes symbols to human gene symbols')
parser.add_argument(
    '--infile',
    type=str,
    help='Input file',
    required=True
)
parser.add_argument(
    '--out',
    type=str,
    help='Output file name',
    required=True
)
args = parser.parse_args()

main(args=args)

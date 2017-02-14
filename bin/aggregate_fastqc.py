"""
This script combines and summarizes the FastQC result files

"""

from crimson import fastqc
import sys
import pandas
import os
import argparse

def main(args):

    data = zip(args.files,args.samples)

    results = {}
    baseq_median = []
    baseq_mean = []
    baseq_10q = []
    baseq_90q = []
    baseq_lower_quartile = []
    baseq_upper_quartile = []
    bases = []

    for d in data:

        assert os.path.splitext(d[0])[1]=='.zip', "Incorrect filetype provided from FastQC"

        os.system('unzip -o -q ' + os.path.abspath(d[0]) + ' -d tmp-fastqc')
        f = fastqc.parse('./tmp-fastqc/stdin_fastqc/fastqc_data.txt')

        results[d[1]] = {
            'Per base sequence quality':f['Per base sequence quality']['status'],
            'Sequence Duplication Levels':f['Sequence Duplication Levels']['status'],
            'Per sequence GC content':f['Per sequence GC content']['status'],
            'Sequence Length Distribution':f['Sequence Length Distribution']['status'],
            'Kmer Content':f['Kmer Content']['status'],
            'Overrepresented sequences':f['Overrepresented sequences']['status'],
            'Per base N content':f['Per base N content']['status'],
            'Per sequence quality scores':f['Per sequence quality scores']['status']
        }
        os.system('rm -rf tmp-fastqc')

        #if args.baseq:

    data = pandas.DataFrame(results)
    data.to_csv(args.out)

parser = argparse.ArgumentParser(description='Aggregates the FastQC results')
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--files',type=str,required=True,nargs='+',help='Space-delimited list of FastQC .zip output files.')
parser.add_argument('--samples',type=str,required=True,nargs='+',help='Space-delimited list of sample names')
parser.add_argument('--baseq',action='store_true',required=False,help='Summarize base quality scores')
args = parser.parse_args()

main(args=args)

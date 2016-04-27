from crimson import fastqc
import sys
import pandas
import os
import argparse

def main(args):

    data = zip(args.files,args.samples)

    results = {}

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
    data = pandas.DataFrame(results)
    data.to_csv(args.out)

parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--files',type=str,required=True,nargs='+')
parser.add_argument('--samples',type=str,required=True,nargs='+')
args = parser.parse_args()

main(args=args)

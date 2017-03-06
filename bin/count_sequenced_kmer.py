import os
import argparse
import sys
from itertools import product

def count_context(args):

    depths = [0,0,0]
    kmer = ['n', 'n', 'n']

    kmers = {''.join(i): 0 for i in product('acgt', repeat=3)}

    for line in sys.stdin:
        if line.find('mpileup') != -1:
            continue
            
        line = line.rstrip().split('\t')
        depth = int(line[3])
        reference = line[2].lower()

        kmer[0] = kmer[1]
        kmer[1] = kmer[2]
        kmer[2] = reference

        depths[0] = depths[1]
        depths[1] = depths[2]
        depths[2] = depth

        if (depths[0] >= args.n) and \
            (depths[1] >= args.n) and \
            (depths[2] >= args.n) and \
            (kmer[0] != 'n') and \
            (kmer[1] != 'n') and \
            (kmer[2] != 'n'):
            kmers[''.join(kmer)] += 1

    fout = open(args.output, 'w')
    fout.write('context,count\n')
    for i in kmers.keys():
        fout.write(i + ',' + str(kmers[i]) + '\n')
    fout.close()

parser = argparse.ArgumentParser(
    description='Count the number of times a certain 3-mer is sequenced'
)
parser.add_argument(
    '--n',
    type=int,
    help='Min depth (default 10)',
    required=True
)
parser.add_argument(
    '--output',
    type=str,
    help='Output file name.',
    required=True
)
args = parser.parse_args()

count_context(args=args)

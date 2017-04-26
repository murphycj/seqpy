import pyfaidx
import pandas
import argparse
import os
import sys
from itertools import islice, product


def window(seq, n=2):
    "Returns a sliding window (of width n) over data from the iterable"
    "   s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   "
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield ''.join(result)
    for elem in it:
        result = result[1:] + (elem,)
        yield ''.join(result)


def main(args):

    if not os.path.exists(args.fasta + '.fai'):
        print '.fai file does not exist!'
        sys.exit()

    fasta = pyfaidx.Fasta(args.fasta)

    kmers = {''.join(i): 0 for i in product('acgt', repeat=3)}

    if args.bed is not None:
        beds = []
        fin = open(args.bed, 'r')
        for line in fin.readlines():
            line = line.rstrip().split('\t')
            assert len(line) >= 3, "found less than 3 entries " + str(line)
            beds.append(line)

        beds = pandas.DataFrame(beds)

        for group in beds.groupby(0):
            chrom = group[0]
            print chrom
            chrom_seq = str(fasta[chrom]).lower()
            for i in group[1].index:
                start = int(group[1].ix[i][1]) - 1
                end = int(group[1].ix[i][2]) - 1
                seq = chrom_seq[start:end]

                for kmer in kmers.keys():
                    kmers[kmer] += seq.count(kmer)

    else:
        for chrom in fasta.keys():
            print chrom
            chrom = str(fasta[chrom]).lower()
            for kmer in kmers.keys():
                kmers[kmer] += chrom.count(kmer)

    fout = open(args.out, 'w')
    for i in kmers.keys():
        fout.write(i + ',' + str(kmers[i]) + '\n')
    fout.close()


parser = argparse.ArgumentParser(
    description='Compute tri-nucleotide frequency from fasta.'
)
parser.add_argument(
    '--fasta',
    type=str,
    required=True,
    help='Fasta file.'
)
parser.add_argument(
    '--kmer',
    type=int,
    required=False,
    help='k-mer length (default 3)'
)
parser.add_argument(
    '--bed',
    type=str,
    required=False,
    default=None,
    help='(Optional) BED file to restrict the analysis to.'
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

main(args=args)

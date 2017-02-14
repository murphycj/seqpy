"""
Convert a MAF file to a matrix with rows as genes and columns as samples.
Each matrix entry just states if the gene is mutant in the sample.

"""

import argparse
import re
import pandas

rr1 = re.compile('^#')
rr2 = re.compile('^Hugo_Symbol')


def main(args):
    fin = open(args.maf, 'r')

    data = {}
    genes = []

    for line in fin.readlines():

        if rr1.findall(line) or rr2.findall(line):
            continue

        line = line.rstrip().split('\t')

        gene = line[0]
        tumor = line[15]

        if gene not in genes:
            genes.append(gene)

        if tumor not in data:
            data[tumor] = [gene]
        else:
            data[tumor].append(gene)

    fin.close()

    for sample, gg in data.items():
        data[sample] = list(set(gg))

    # write output

    samples = data.keys()
    samples.sort()
    genes.sort()

    fout = open(args.out, 'w')
    fout.write('gene,' + ','.join(samples) + '\n')

    for g in genes:
        fout.write(g)
        for s in samples:
            if g in data[s]:
                fout.write(',MUT')
            else:
                fout.write(',-')
        fout.write('\n')
    fout.close()


parser = argparse.ArgumentParser(description='MAF to matrix')

parser.add_argument(
    '--maf',
    type=str,
    required=True,
    help='MAF file.'
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

main(args=args)

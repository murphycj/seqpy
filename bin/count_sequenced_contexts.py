"""

"""

import pandas
import os
import argparse
from subprocess import Popen, PIPE
import vcf


SUBSTITUTIONS = {
    'a>c':[],
    'a>g':[],
    'a>t':[],
    'g>a':[],
    'g>c':[],
    'g>t':[]
}

CONTEXTS = []
for i in ['a','c','t','g']:
    for j in ['a','c','t','g']:
        CONTEXTS.append([i,j])

#12 possible substitutions, which can be effectively reduced to 6
CONVERSIONS = {
    'a':'a',
    'g':'g',
    't':'a',
    'c':'g'
}


def count_context(args):

    depth = 0
    contexts = ['','','']
    three = 0
    bases = 0

    context_counts = []
    for i in CONTEXTS:
        for j in SUBSTITUTIONS.keys():
            context_counts.append(i[0] + j[0] + i[1])
    context_counts = {i:0 for i in context_counts}

    p = Popen([args.samtools,'mpileup','-f',args.genome,args.bam],stdout=PIPE,stderr=PIPE,stdin=PIPE)

    for line in p.stdout.readlines():
        line = line.rstrip().split('\t')
        depth = int(line[3])
        reference = line[2].lower()

        contexts[0] = contexts[1]
        contexts[1] = contexts[2]

        if depth>=1 and reference!='n':
                three += 1
                contexts[2] = reference
        else:
            three = 0
            contexts[2] = ''

        if three > 3:
            three = 3

        if three == 3:
            context = contexts[0] + CONVERSIONS[contexts[1]] + contexts[2]
            context_counts[context] += 1

    fout = open(args.prefix + '.context_count.txt','w')
    fout.write('context,count\n')
    for i in context_counts.keys():
        fout.write(i + ',' + str(context_counts[i]) + '\n')
    fout.close()

    return context_counts

def main(args):

    assert os.path.isfile(args.genome), "Genome does not exist"
    assert os.path.isfile(args.genome + '.fai'), "Index your genome"

    context_counts = count_context(args=args)


    data.to_csv(args.out)


parser = argparse.ArgumentParser(description='Get the mutation context spectrum')
parser.add_argument('--bam',type=str,help='BAM file.',required=True)
parser.add_argument('--genome',type=str,help='Path to genome (e.g. /path/to/mm9.fa).',required=True)
parser.add_argument('--samtools',type=str,help='Path to samtools',required=True)
parser.add_argument('--output',type=str,help='Output file name.',required=True)
args = parser.parse_args()

main(args=args)

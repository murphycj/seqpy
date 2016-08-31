import os
import argparse
import sys


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

    for line in sys.stdin:
        if line.find('mpileup')!=-1:
            continue
        try:
            line = line.rstrip().split('\t')
            depth = int(line[3])
            reference = line[2].lower()

            contexts[0] = contexts[1]
            contexts[1] = contexts[2]

            if depth>=args.n and reference!='n':
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
        except:
            import pdb; pdb.set_trace()

    fout = open(args.output,'w')
    fout.write('context,count\n')
    for i in context_counts.keys():
        fout.write(i + ',' + str(context_counts[i]) + '\n')
    fout.close()

parser = argparse.ArgumentParser(description='Count the number of times a certain tri-nucleotide is sequenced')
parser.add_argument('--n',type=int,help='Min depth (default 10)',required=True)
parser.add_argument('--output',type=str,help='Output file name.',required=True)
args = parser.parse_args()

count_context(args=args)

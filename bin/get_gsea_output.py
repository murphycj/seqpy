import os
import sys
import argparse
from seqpy import parsers
import glob


def main(args):

    if args.zip:
        os.mkdir('tmpGSEA')
        os.system('unzip ' + args.indir + ' -d tmpGSEA')
        g = parsers.GSEA('./tmpGSEA', args.out, args.phenotype1, args.phenotype2)
    elif args.targz:
        os.mkdir('tmpGSEA')
        os.system('tar -xf ' + args.indir + ' -C tmpGSEA')
        directory = glob.glob('./tmpGSEA/*')
        assert len(directory)==1,"did not find 1 directory!"
        g = parsers.GSEA(directory[0], args.out, args.phenotype1, args.phenotype2)
    else:
        g = parsers.GSEA(args.indir, args.out, args.phenotype1, args.phenotype2)


    if args.phenotype_1_tab_name is not None:
        phenotype1 = args.phenotype_1_tab_name
    else:
        phenotype1 = args.phenotype1

    if args.phenotype_2_tab_name is not None:
        phenotype2 = args.phenotype_2_tab_name
    else:
        phenotype2 = args.phenotype2

    g.parse_pathway_excel(
        leadingEdge=args.leadingEdge,
        leadingEdgeGenes=args.leadingEdgeGenes,
        allGenes=args.allGenes
    )

    g.write_to_excel(args.out + '.xlsx',phenotype1=phenotype1,phenotype2=phenotype2)

    if args.zip:
        os.system('rm -rf tmpGSEA')

parser = argparse.ArgumentParser(description='Get the gsea output.')
parser.add_argument(
    '--indir',
    type=str,
    help='directory containing results',
    required=True
)
parser.add_argument(
    '--zip',
    action='store_true',
    help='(Optional) The provided GSEA output directory is a zip file ',
    required=False
)
parser.add_argument(
    '--targz',
    action='store_true',
    help='(Optional) The provided GSEA output directory is a tar.gz file ',
    required=False
)
parser.add_argument(
    '--phenotype1',
    type=str,
    required=False,
    default='up',
    help='Phenotype 1 (default up)'
)
parser.add_argument(
    '--phenotype2',
    type=str,
    required=False,
    default='down',
    help='Phenotype 2 (default down)'
)
parser.add_argument(
    '--phenotype_1_tab_name',
    type=str,
    default=None,
    help='(Optional) Tab name for first phenotype results.',
    required=False
)
parser.add_argument(
    '--phenotype_2_tab_name',
    type=str,
    default=None,
    help='(Optional) Tab name for second phenotype results.',
    required=False
)
parser.add_argument(
    '--leadingEdge',
    action='store_true',
    help='(Optional) Include in the output the number of genes in the ' +
         'leading edge gene set.',
    required=False
)
parser.add_argument(
    '--leadingEdgeGenes',
    action='store_true',
    help='(Optional) Include in the output all the genes in the leading ' +
         'leading edge gene set.',
    required=False
)
parser.add_argument(
    '--allGenes',
    action='store_true',
    help='(Optional) Include in the output the genes in the ' +
         ' gene set.',
    required=False
)
parser.add_argument('--out', type=str, help='Output prefix', required=True)
args = parser.parse_args()

main(args=args)

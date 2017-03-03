import os
import sys
import argparse
from seqpy import parsers


def main(args):

    if args.zip:
        os.mkdir('tmpGSEA')
        os.system('unzip ' + args.indir + ' -d tmpGSEA')
        g = parsers.GSEA('./tmpGSEA', args.out)
    else:
        g = parsers.GSEA(args.indir, args.out)

    if args.reverse:
        g.parse_pathway_excel(
            reverse=True,
            leadingEdge=args.leadingEdge,
            leadingEdgeGenes=args.leadingEdgeGenes,
            allGenes=args.allGenes
        )
    else:
        g.parse_pathway_excel(
            reverse=False,
            leadingEdge=args.leadingEdge,
            leadingEdgeGenes=args.leadingEdgeGenes,
            allGenes=args.allGenes
        )

    g.write_to_excel(args.out + '.xlsx')

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
    '--reverse',
    action='store_true',
    help='(Optional) By default the gsea_report_for_na_neg_* file ' +
         'regarded at up-regulate genes',
    required=False
)
parser.add_argument(
    '--up_tab_name',
    type=str,
    help='(Option) Tab name for up-regualted genes.',
    required=False
)
parser.add_argument(
    '--down_tab_name',
    type=str,
    help='(Optional) Tab name for down-regualted genes.',
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

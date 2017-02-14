import argparse
from seqpy import gsea


def main(args):
    g = gsea.GSEA(args.indir,args.out)

    if args.reverse:
        g.parse_pathway_excel(reverse=True)
    else:
        g.parse_pathway_excel()
    g.write_to_excel(args.out + '.xlsx')


parser = argparse.ArgumentParser(description='Get the gsea output.')
parser.add_argument('--indir',type=str,help='directory containing results',required=True)
parser.add_argument('--reverse',action='store_true',help='By default the gsea_report_for_na_neg_* file regarded at up-regulate gebes',required=False)
parser.add_argument('--up_tab_name',type=str,help='Tab name for up-regualted genes.',required=False)
parser.add_argument('--down_tab_name',type=str,help='Tab name for down-regualted genes.',required=False)
parser.add_argument('--out',type=str,help='Output prefix',required=True)
args = parser.parse_args()

main(args=args)

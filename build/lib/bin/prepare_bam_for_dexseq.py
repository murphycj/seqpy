"""
Prepares the input BAM file for DEXseq script processing

"""

import argparse
import pysam


def main(args):

    if pysam.__version__ == '0.7.5':
        bamfile = pysam.Samfile(args.bam,'rb')
        out_bamfile = pysam.Samfile(args.out,'wb',template=bamfile)
    else:
        bamfile = pysam.AlignmentFile(args.bam,'rb')
        out_bamfile = pysam.AlignmentFile(args.out,'wb',template=bamfile)
    tag = ('XT','U')
    tag_add = [('NH',1)]

    for alignment in bamfile:
        if args.XTtag:
            if alignment.mapq >= args.mapq and tag in alignment.tags:
                alignment.tags = alignment.tags + tag_add
                out_bamfile.write(alignment)
        else:
            if alignment.mapq >= args.mapq: # and tag in alignment.tags:
                alignment.tags = alignment.tags + tag_add
                out_bamfile.write(alignment)

    out_bamfile.close()
    bamfile.close()


parser = argparse.ArgumentParser(description='Prepares the input BAM file for DEXseq script processing')
parser.add_argument('--bam',type=str,help='BAM file',required=True)
parser.add_argument('--out',type=str,help='output BAM file name',required=True)
parser.add_argument('--mapq',type=int,help='MAPQ minimum (default 20)',required=False,default=20)
parser.add_argument('--XTtag',help='Include if alignments need the XT:i:U tag.',action='store_true')

args = parser.parse_args()

main(args=args)

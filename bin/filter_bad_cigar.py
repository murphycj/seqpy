import pysam
import argparse

def main(args):

    samfile = pysam.AlignmentFile(args.infile, "rb")
    outfile = pysam.AlignmentFile(args.out, "wb", template=samfile)

    n = 0

    for read in samfile:
        if read.cigarstring.find('M')==-1:
            n+=1
        else:
            outfile.write(read)

    samfile.close()
    outfile.close()

    print "Found %s reads without M!" % n


parser = argparse.ArgumentParser(description='Filters mapped reads that have no M in the CIGAR string')
parser.add_argument(
    '--infile',
    type=str,
    required=True,
    help='Input BAM file'
)

parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output BAM file'
)
args = parser.parse_args()

main(args=args)

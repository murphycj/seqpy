import vcf
import sys
import os
import argparse

def main(args):

    if not os.path.exists(args.sample):
        os.mkdir(args.sample)

    vcf_in = vcf.Reader(open(args.vcf,'r'))
    fout = open(args.sample + '.igv_batch.sh','w')
    fout.write('new\n')
    fout.write('genome ' + args.genome + '\n')
    fout.write('load ' + args.bam + '\n')
    fout.write('snapshotDirectory .\n')

    for v in vcf_in:

        for s in v.samples:
            if s.sample==args.sample and s.called:
                snapshot = './' + args.sample + '/' + args.sample + '.chr' + str(v.CHROM) + '.' + str(v.POS) + '.png'

                fout.write('goto chr' + str(v.CHROM) + ':' + str(v.POS) + '\n')

                if args.collapse:
                    fout.write('collapse\n')

                fout.write('snapshot ' + snapshot + '\n')

    fout.write('exit\n')
    fout.close()

parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,help='VCF file',required=True)
parser.add_argument('--sample',type=str,help='samples',required=True)
parser.add_argument('--bam',type=str,help='BAM',required=True)
parser.add_argument('--genome',type=str,help='The genome to load (e.g. hg19, mm9, mm10, ...)',required=True)
parser.add_argument('--collapse',action='store_true',help='Whether to collapse the reads.',required=False)
args = parser.parse_args()

main(args=args)

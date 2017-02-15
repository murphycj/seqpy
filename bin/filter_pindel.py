"""
Filter a VCF file output from pindel tumor-normal mutation calling

"""

import vcf
import argparse


def main(args):
    vcf_in = vcf.Reader(open(args.vcf, 'r'))
    vcf_out = vcf.Writer(open(args.out, 'w'), vcf_in)

    for v in vcf_in:
        keep = False

        # get sample index

        samples = map(lambda x: x.sample, v.samples)
        n = 0
        for i in samples:
            if i == args.tumor:
                tumor_index = n
            if i == args.normal:
                normal_index = n
            n += 1

        # check minimum coverage in normal

        if hasattr(v.samples[normal_index].data, 'AD'):
            AD = v.samples[normal_index].data.AD
            if sum(AD) != 0:
                vaf = AD[1] / float(sum(AD))
            else:
                vaf = 0.0
            if sum(AD) < args.min_normal or vaf > args.max_normal_vaf:
                continue
            if args.max_normal_AD != -1 and AD[1] > args.max_normal_AD:
                continue
        else:
            print 'Site does not have AD attr'
            print v

        # check minimum coverage in tumor

        if hasattr(v.samples[tumor_index].data, 'AD'):
            AD = v.samples[tumor_index].data.AD
            if sum(AD) < args.min_tumor:
                continue
            freq = AD[1] / float(sum(AD))
            if AD[1] < args.AD or freq < args.freq:
                continue
        else:
            print 'Site does not have AD attr'
            print v

        vcf_out.write_record(v)
    vcf_out.close()


parser = argparse.ArgumentParser(description='Filter Pindel output')
parser.add_argument(
    '--vcf',
    type=str,
    help='VCf file',
    required=True
)
parser.add_argument(
    '--tumor',
    type=str,
    help='Tumor sample name (default TUMOR)',
    default='TUMOR',
    required=False
)
parser.add_argument(
    '--normal',
    type=str,
    help='Tumor sample name (default NORMAL)',
    default='NORMAL',
    required=False
)
parser.add_argument(
    '--AD',
    type=int,
    help='Min mutant-supporting reads in tumor (default 4)',
    default=4,
    required=False
)
parser.add_argument(
    '--freq',
    type=float,
    help='Min frequency (default 0.05)',
    default=0.05,
    required=True
)
parser.add_argument(
    '--min_tumor',
    type=int,
    help='Min coverage in tumor (default 12)',
    default=12,
    required=False
)
parser.add_argument(
    '--min_normal',
    type=int,
    help='Min coverage in normal (default 12)',
    default=12,
    required=False
)
parser.add_argument(
    '--max_normal_vaf',
    type=float,
    help='Max variant allele frequency in normal (default 0.01)',
    default=0.01,
    required=False
)
parser.add_argument(
    '--max_normal_AD',
    type=int,
    help='Maximum number of mutant-supporting reads in normal (default none)',
    default=-1,
    required=False
)
parser.add_argument(
    '--out',
    type=str,
    help='out vcf file',
    required=True
)
args = parser.parse_args()

main(args=args)

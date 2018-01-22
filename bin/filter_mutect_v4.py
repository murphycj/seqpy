"""
Filter a VCF file output from mutect tumor-normal mutation calling

"""

import vcf
import argparse


def main(args):
    vcf_in = vcf.Reader(open(args.vcf, 'r'))
    vcf_out = vcf.Writer(open(args.out, 'w'), vcf_in)

    for v in vcf_in:

        TLOD = v.INFO['TLOD']
        keep = False
        if type(TLOD)==list:
            for i in TLOD:
                if float(i) >= args.tlod:
                    keep=True
        else:
            if float(TLOD) >= args.tlod:
                keep=True

        if not keep:
            continue

        if not args.nonpaired:
            NLOD = v.INFO['NLOD']
            keep = False
            if type(NLOD)==list:
                for i in NLOD:
                    if float(i)>=args.nlod:
                        keep=True
            else:
                if float(NLOD) >= args.nlod:
                    keep=True

            if not keep:
                continue

        # get sample index

        samples = map(lambda x: x.sample, v.samples)
        n = 0
        for i in samples:
            if i == args.tumor:
                tumor_index = n
            if not args.nonpaired and i == args.normal:
                normal_index = n
            n += 1

        # check minimum coverage in normal

        if not args.nonpaired:
            if hasattr(v.samples[normal_index].data, 'AD'):
                AD = v.samples[normal_index].data.AD
                if sum(AD) < args.min_normal:
                    continue
                if AD[1] >= args.max_normal_support:
                    continue
            else:
                print 'Site does not have AD attr'
                print v

        # check minimum coverage in tumor

        if hasattr(v.samples[tumor_index].data, 'AD'):
            AD = v.samples[tumor_index].data.AD
            freq = float(v.samples[tumor_index].data.AF)
            if AD[1] < args.AD or sum(AD) < args.min_tumor or freq < args.freq:
                continue
        else:
            print 'Site does not have AD attr'
            print v

        vcf_out.write_record(v)
    vcf_out.close()


parser = argparse.ArgumentParser(description='')
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
    '--nonpaired',
    action='store_true',
    help='Mutect file is from running on single sample (not tumor-normal paired).',
    required=False
)
parser.add_argument(
    '--AD',
    type=int,
    help='Min alt-supporting reads (default 4)',
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
    '--max_normal_support',
    type=int,
    help='Max alt allele support in normal (default 3).',
    default=3,
    required=False
)
parser.add_argument(
    '--tlod',
    type=int,
    help='Min TLOD',
    default=6.3,
    required=False
)
parser.add_argument(
    '--nlod',
    type=int,
    help='Min NLOD',
    default=2.2,
    required=False
)
parser.add_argument(
    '--artlod',
    type=int,
    help='Max artifact LOD',
    default=2.2,
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

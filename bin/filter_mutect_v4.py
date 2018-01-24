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

        if len(TLOD) > 2:
            print('multiple TLOD!!')
            print(v)
            exit()

        if type(TLOD) == list:
            TLOD = TLOD[0]

        if float(TLOD) >= args.tlod:
            continue

        if not args.nonpaired:
            NLOD = v.INFO['NLOD']

            if len(NLOD) > 2:
                print('multiple NLOD!!')
                print(v)
                exit()

            if type(NLOD) == list:
                NLOD = NLOD[0]

            if float(NLOD) < args.nlod:
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

        # check coverage in normal

        if not args.nonpaired:
            if hasattr(v.samples[normal_index].data, 'AD'):
                AD = v.samples[normal_index].data.AD
                coverage = float(sum(AD))

                if coverage < args.min_normal:
                    continue

                if len(AD) > 2:
                    print('multiple AD!!')
                    print(v)
                    exit()
                else:
                    if args.max_normal_support!=-1 and AD[1] >= args.max_normal_support:
                        continue
                    if float(AD[1])/coverage >= args.max_normal_af:
                        continue
            else:
                print('Site does not have AD attr')
                print(v)
                exit()

        # check coverage in tumor

        if hasattr(v.samples[tumor_index].data, 'AD'):
            AD = v.samples[tumor_index].data.AD
            coverage = float(sum(AD))

            if coverage < args.min_tumor:
                continue

            if len(AD) > 2:
                print('multiple AD!!')
                print(v)
                exit()
            else:
                if AD[1] < args.AD or AD[1]/coverage < args.freq:
                    continue
        else:
            print('Site does not have AD attr')
            print(v)
            exit()

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
    help='Max alt allele support in normal (default -1).',
    default=-1,
    required=False
)
parser.add_argument(
    '--max_normal_af',
    type=float,
    help='Max alt allele AF in normal (default 0.03).',
    default=0.03,
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

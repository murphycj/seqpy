"""
Filter a VCF file output from mutect tumor-normal mutation calling

"""

import vcf
import argparse


def main(args):
    vcf_in = vcf.Reader(open(args.vcf, 'r'))
    vcf_out = vcf.Writer(open(args.out, 'w'), vcf_in)

    for v in vcf_in:
        keep = False

        if args.homologous_mapping_event and 'homologous_mapping_event' in v.FILTER:
            continue

        if args.clustered_events and 'clustered_events' in v.FILTER:
            continue

        if args.triallelic_site and 'triallelic_site' in v.FILTER:
            continue

        if args.alt_allele_in_normal and 'alt_allele_in_normal' in v.FILTER:
            continue

        if args.multi_event_alt_allele_in_normal and 'multi_event_alt_allele_in_normal' in v.FILTER:
            continue

        if float(v.INFO['TLOD']) < args.tlod:
            continue

        if float(v.INFO['NLOD']) < args.nlod:
            continue

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
            if sum(AD) < args.min_normal:
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
    '--tlod',
    type=int,
    help='Min TLOD',
    default=6.3,
    required=False
)
parser.add_argument(
    '--nlod',
    type=int,
    help='Min TLOD',
    default=2.2,
    required=False
)
parser.add_argument(
    '--out',
    type=str,
    help='out vcf file',
    required=True
)
parser.add_argument(
    '--homologous_mapping_event',
    action='store_true',
    help='Remove homologous_mapping_event',
    required=False
)
parser.add_argument(
    '--clustered_events',
    action='store_true',
    help='Remove clustered_events',
    required=False
)
parser.add_argument(
    '--triallelic_site',
    action='store_true',
    help='Remove triallelic_site',
    required=False
)
parser.add_argument(
    '--alt_allele_in_normal',
    action='store_true',
    help='Remove alt_allele_in_normal',
    required=False
)
parser.add_argument(
    '--multi_event_alt_allele_in_normal',
    action='store_true',
    help='Remove multi_event_alt_allele_in_normal',
    required=False
)
args = parser.parse_args()

main(args=args)

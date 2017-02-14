"""
Takes a Mutect vcf from a single sample and removes one of the
QSS field so bcftools merge can be used
"""
import argparse
import re
import vcf


def main(args):

    vcf_in = vcf.Reader(open(args.vcf, 'r'))
    vcf_out = vcf.Writer(open(args.out, 'w'),vcf_in)

    for variant in vcf_in:
        for sample in variant:
            sample['QSS'].pop()
        vcf_out.write_record(variant)


parser = argparse.ArgumentParser(
    description='Takes a Mutect vcf from a single sample ' +
    'and removes one of the QSS field so bcftools merge can be used'
)
parser.add_argument(
    '--vcf',
    type=str,
    help='VCF file',
    required=True
)
parser.add_argument(
    '--out',
    type=str,
    help='out results',
    required=False
)
args = parser.parse_args()

main(args=args)

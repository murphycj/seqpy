import vcf
import argparse

def main(args):
    vcf_in = vcf.Reader(open(args.vcf,'r'))
    vcf_out = vcf.Writer(open(args.out,'w'), vcf_in)

    for v in vcf_in:
        keep = False
        for s in v.samples:
            if s.called and s.data.AD>=args.alt:
                keep=True
                break

        if keep:
            vcf_out.write_record(v)
    vcf_out.close()

parser = argparse.ArgumentParser(description='Given a VCf file of control samples, keep only those sites with at least one conrol with n alt-supporting reads')
parser.add_argument('--vcf',type=str,help='VCf file',required=True)
parser.add_argument('--alt',type=int,help='Min alt-supporting reads',required=True)
parser.add_argument('--out',type=str,help='out vcf file',required=False)
args = parser.parse_args()

main(args=args)

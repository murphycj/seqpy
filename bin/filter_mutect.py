import vcf
import argparse

def main(args):
    vcf_in = vcf.Reader(open(args.vcf,'r'))
    vcf_out = vcf.Writer(open(args.out,'w'), vcf_in)

    for v in vcf_in:
        keep = False

        if args.homologous_mapping_event and 'homologous_mapping_event' in v.FILTER:
            continue

        if args.clustered_events and 'clustered_events' in v.FILTER:
            continue

        if args.triallelic_site and 'triallelic_site' in v.FILTER:
            continue

        if float(v.INFO['TLOD']) < args.tlod:
            continue

        if float(v.INFO['NLOD']) < args.nlod:
            continue


        if hasattr(v.samples[0].data,'AD'):
            AD = v.samples[0].data.AD
            freq = float(v.samples[0].data.AF)
            if AD[1] < args.AD or sum(AD) < args.min or freq < args.freq:
                continue
        else:
            continue

        vcf_out.write_record(v)
    vcf_out.close()

parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,help='VCf file',required=True)
parser.add_argument('--AD',type=int,help='Min alt-supporting reads',required=True)
parser.add_argument('--freq',type=float,help='Min frequency',required=True)
parser.add_argument('--min',type=int,help='Min reads',required=True)
parser.add_argument('--tlod',type=int,help='Min TLOD',required=True)
parser.add_argument('--nlod',type=int,help='Min TLOD',required=True)
parser.add_argument('--out',type=str,help='out vcf file',required=False)
parser.add_argument('--homologous_mapping_event',action='store_true',help='Remove homologous_mapping_event',required=False)
parser.add_argument('--clustered_events',action='store_true',help='Remove clustered_events',required=False)
parser.add_argument('--triallelic_site',action='store_true',help='Remove triallelic_site',required=False)
args = parser.parse_args()

main(args=args)

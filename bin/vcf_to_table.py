"""
Convert a vcf file annotated by SnpEff into a csv file that has one mutation
per line (or multiple lines for that particular mutation in case its effect
has varying context) and with samples as columns to indicate which samples
have which mutations
"""

import vcf
import argparse

def main(args):
    vcf_in = vcf.Reader(open(args.vcf,'r'))
    fout = open(args.out,'w')
    fout.write('CHROM,POS,REF,ALT,GENE,TRANSCRIPT,BASE PAIR CHANGE,AMINO ACID CHANGE,EFFECT,SAMPLE COUNT')
    for s in vcf_in.samples:
        fout.write(',' + s)
    fout.write('\n')

    for v in vcf_in:

        #loop through all possible mutation annotations for each mutation

        genes = []
        if 'ANN' not in v.INFO:
            continue

        #get the total number of samples with the mutation

        n=0
        for s in v.samples:
            if s.called:
                n+=1

        for ann in v.INFO['ANN']:

            #split the ANN fields and get predicted effect

            info = ann.split('|')
            if (not args.everything) and (info[1] not in args.filter_effects):
                continue

            #if has effect that is desired, print it

            fout.write(
                str(v.CHROM) + ',' +
                str(v.POS) + ',' +
                v.REF + ',' +
                str(v.ALT[0]) + ',' +
                str(info[3]) + ',' +
                info[6] + ',' +
                info[9] + ',' +
                info[10] + ',' +
                info[1]
            )

            fout.write(',' + str(n))

            for s in v.samples:

                #if mutation is present, print the mutation frequency

                if s.called:
                    fout.write(',%s (%s/%s)' % (s.data.FREQ, s.data.AD, (s.data.AD + s.data.RD)))
                    n+=1
                else:
                    fout.write(',-')
            fout.write('\n')
    fout.close()


parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,help='Meta-directory where each sub-directory contains fpkm for each sample',required=True)
parser.add_argument(
    '--filter_effects',
    type=str,
    nargs='+',
    required=False,
    help='Space-delimited list of variant types to include in output (default: ' + \
    'missense_variant frameshift_variant start_lost stop_gained ' + \
    '5_prime_UTR_premature_start_codon_gain_variant nonsense_mediated_decay ' + \
    'inframe_insertion disruptive_inframe_insertion inframe_deletion disruptive_inframe_deletion ' + \
    'rare_amino_acid_variant splice_acceptor_variant splice_donor_variant ' + \
    'stop_lost non_coding_exon_variant)',
    default = [
        'missense_variant',
        'frameshift_variant',
        'start_lost',
        'stop_gained',
        '5_prime_UTR_premature_start_codon_gain_variant',
        'nonsense_mediated_decay',
        'inframe_insertion',
        'disruptive_inframe_insertion',
        'inframe_deletion',
        'disruptive_inframe_deletion',
        'rare_amino_acid_variant',
        'splice_acceptor_variant',
        'splice_donor_variant',
        'stop_lost',
        'non_coding_exon_variant'
        ]
)
parser.add_argument('--everything',action='store_true',help='Output all types of variants',required=False)
parser.add_argument('--out',type=str,help='Collate FPKM values for genes',required=False)
args = parser.parse_args()

main(args=args)

"""
Convert a vcf file annotated by SnpEff into a csv file that has one mutation
per line (or multiple lines for that particular mutation in case its effect
has varying context) and with samples as columns to indicate which samples
have which mutations
"""

import vcf
import argparse
from seqpy.parsers import SnpEff, SnpEffInfo


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

        vinfo = SnpEffInfo(v.INFO)

        if not vinfo.has_ann():
            continue

        #get the total number of samples with the mutation

        n=0
        for s in v.samples:
            if s.called:
                n+=1

        for ann in vinfo.ann:

            if (not args.everything):
                has_right_annotation=False

                for a in ann.annotation:
                    if a in args.filter_effects:
                        has_right_annotation=True

                if not has_right_annotation:
                    continue

            #if has effect that is desired, print it

            fout.write(
                str(v.CHROM) + ',' +
                str(v.POS) + ',' +
                v.REF + ',' +
                str(v.ALT[0]) + ',' +
                str(ann.gene_name) + ',' +
                ann.feature_id + ',' +
                ann.basepair_change + ',' +
                ann.aminoacid_change + ',' +
                '&'.join(ann.annotation)
            )

            fout.write(',' + str(n))

            for s in v.samples:

                #if mutation is present, print the mutation frequency

                if s.called:
                    if args.mutect:
                        fout.write(',%s (%s/%s)' % (s.data.AF, s.data.AD[1], sum(s.data.AD)))
                    else:
                        fout.write(',%s (%s/%s)' % (s.data.FREQ, s.data.AD, (s.data.AD + s.data.RD)))
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
    'stop_lost non_coding_exon_variant). Use --filter_effects none if you wanted to include all effects',
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
parser.add_argument('--mutect',action='store_true',help='input vcf is from Mutect',required=False)
parser.add_argument('--out',type=str,help='Collate FPKM values for genes',required=False)
args = parser.parse_args()

main(args=args)

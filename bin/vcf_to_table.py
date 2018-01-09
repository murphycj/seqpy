"""
Convert a vcf file annotated by SnpEff into a csv file that has one mutation
per line (or multiple lines for that particular mutation in case its effect
has varying context) and with samples as columns to indicate which samples
have which mutations
"""

import re
import argparse
from collections import Counter

import vcf
from seqpy.parsers.SNPEff import SnpEffInfo
from seqpy.parsers.VCF import VAF


def main(args):

    cc = re.compile('^ENS')

    vcf_in = vcf.Reader(open(args.vcf, 'r'))
    fout = open(args.out, 'w')

    if args.cosmic:
        fout.write('CHROM,POS,REF,ALT,GENE,TRANSCRIPT,BASE PAIR CHANGE,AMINO ACID CHANGE,EFFECT,COSMIC,SAMPLE COUNT,FLAGS')
    else:
        fout.write('CHROM,POS,REF,ALT,GENE,TRANSCRIPT,BASE PAIR CHANGE,AMINO ACID CHANGE,EFFECT,SAMPLE COUNT,FLAGS')

    for sample in vcf_in.samples:
        fout.write(',' + sample)
    fout.write('\n')

    for v in vcf_in:

        # loop through all possible mutation annotations for each mutation

        genes = []

        vinfo = SnpEffInfo(v.INFO)

        if args.cosmic:
            if v.ID is None:
                cosmic = '-'
            else:
                cosmic = str(v.ID)
        else:
            cosmic = '-'

        if not vinfo.has_ann():
            print 'No annotation, skipping: ' + str(v)
            continue

        # get the total number of samples with the mutation

        n = 0
        for sample in v.samples:
            if sample.called:
                n += 1

        # handles the case that there are multiple annotations for a
        # gene even though using only canonical

        annotations = {}
        for ann in vinfo.ann:
            if ann.gene_id not in annotations:
                annotations[ann.gene_id] = '&'.join(ann.annotation)
            else:
                annotations[ann.gene_id] += '&' + '&'.join(ann.annotation)

        for ann in vinfo.ann:

            if (not args.everything):
                has_right_annotation = False

                for a in ann.annotation:
                    if a in args.filter_effects:
                        has_right_annotation = True

                if not has_right_annotation:
                    continue

            if not args.NET:
                if not cc.findall(ann.feature_id):
                    continue

            if ann.gene_id not in annotations:
                continue

            # if has effect that is desired, print it

            fout.write(
                str(v.CHROM) + ',' +
                str(v.POS) + ',' +
                v.REF + ',' +
                str(v.ALT[0]) + ',' +
                str(ann.gene_name) + ',' +
                ann.feature_id + ',' +
                ann.basepair_change + ',' +
                ann.aminoacid_change + ',' +
                annotations[ann.gene_id] + ',' +
                ';'.join(v.FILTER)
            )

            del annotations[ann.gene_id]

            if args.cosmic:
                fout.write(',' + cosmic + ',' + str(n))
            else:
                fout.write(',' + str(n))

            for sample in v.samples:

                # if mutation is present, print the mutation frequency

                if sample.called:
                    vaf = VAF(
                        sample=sample,
                        mutect=args.mutect,
                        varscan=args.varscan,
                        pindel=args.pindel
                    )
                    try:
                        fout.write(
                            ',%s (%s/%s)' %
                            (
                                round(vaf.freq, 4),
                                vaf.mutant,
                                vaf.mutant + vaf.reference,
                            )
                        )
                    except:
                        import pdb; pdb.set_trace()
                else:
                    fout.write(',-')
            fout.write('\n')
    fout.close()


parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--vcf',
    type=str,
    help='Meta-directory where each sub-directory contains fpkm for ' +
         'each sample',
    required=True
)
parser.add_argument(
    '--filter_effects',
    type=str,
    nargs='+',
    required=False,
    help='Space-delimited list of variant types to include in output ' +
    '(default: missense_variant frameshift_variant start_lost stop_gained ' +
    '5_prime_UTR_premature_start_codon_gain_variant nonsense_mediated_decay ' +
    'inframe_insertion disruptive_inframe_insertion inframe_deletion ' +
    'disruptive_inframe_deletion rare_amino_acid_variant ' +
    'splice_acceptor_variant splice_donor_variant ' +
    'stop_lost non_coding_exon_variant). Use --filter_effects ' +
    'none if you wanted to include all effects',
    default=[
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
parser.add_argument(
    '--everything',
    action='store_true',
    help='Output all types of variants',
    required=False
)
parser.add_argument(
    '--mutect',
    action='store_true',
    help='Input vcf is from Mutect (this and/or --varscan must be set)',
    required=False
)
parser.add_argument(
    '--varscan',
    action='store_true',
    help='Input vcf is from varscan (this and/or --mutect must be set)',
    required=False
)
parser.add_argument(
    '--pindel',
    action='store_true',
    help='Input vcf contains indel calls from Pindel.',
    required=False
)
parser.add_argument(
    '--cosmic',
    action='store_true',
    help='ID column contains COSMIC mutations',
    required=False
)
parser.add_argument(
    '--NET',
    action='store_true',
    help='Include non-ensembl transcripts',
    required=False
)
parser.add_argument(
    '--out',
    type=str,
    help='Out file name',
    required=False
)
args = parser.parse_args()

main(args=args)

import vcf
import argparse


EFFECTs = [
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

        for ann in v.INFO['ANN']:
            info = ann.split('|')
            if info[1] not in EFFECTs:
                continue


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
            n=0
            for s in v.samples:
                if s.called:
                    n+=1
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
parser.add_argument('--out',type=str,help='Collate FPKM values for genes',required=False)
args = parser.parse_args()

main(args=args)

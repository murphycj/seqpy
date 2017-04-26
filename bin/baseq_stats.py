import sys
import numpy as np
import pandas
import os
import argparse
from Bio import SeqIO
import ggplot

def main(args):

    if args.fastq is not None:
        sample = args.sample
        q_scores = list()

        for fastq in args.fastq:
            if args.gz:
                os.system('gunzip -c ' + fastq + " > basq_stats.tmp.fastq")
                fastq_file = 'basq_stats.tmp.fastq'
            else:
                fastq_file = fastq

            seqio = SeqIO.parse(fastq_file, format='fastq')

            for seq_record in seqio:
                q_scores += seq_record.letter_annotations['phred_quality']

        df = pandas.DataFrame(q_scores,columns=['BaseQuality'])

        p = ggplot.ggplot(ggplot.aes(x='BaseQuality'), data=df) + ggplot.geom_histogram(binwidth=1) + ggplot.labs(title = sample)
        p.save(args.prefix + '-hist.png')

        # statistics
        fout = open(args.prefix + '-statistics.csv','w')
        for i in np.arange(0.1,1,0.1):
            fout.write(str(int(100*i)) + '% quantile,' + str(df.quantile(i)[0]) + '\n')
        for i in np.arange(1,41):
            fout.write('Proportion bases >= Q' + str(i) + ',' + str((df>=i).sum()[0]/float(df.shape[0])) + '\n')

        fout.close()

        if os.path.exists('basq_stats.tmp.fastq'):
            os.system('rm basq_stats.tmp.fastq')
    else:
        import pysam
        sample = args.sample
        q_scores = list()

        bam = pysam.AlignmentFile(args.bam,'rb')

        for seq_record in bam:
            q_scores += pysam.qualitystring_to_array(seq_record.qual).tolist()

        df = pandas.DataFrame(q_scores,columns=['BaseQuality'])

        p = ggplot.ggplot(ggplot.aes(x='BaseQuality'), data=df) + ggplot.geom_histogram(binwidth=1) + ggplot.labs(title = sample)
        p.save(args.prefix + '-hist.png')

        # statistics
        fout = open(args.prefix + '-statistics.csv','w')
        for i in np.arange(0.1,1,0.1):
            fout.write(str(int(100*i)) + '% quantile,' + str(df.quantile(i)[0]) + '\n')
        for i in np.arange(1,41):
            fout.write('Proportion bases >= Q' + str(i) + ',' + str((df>=i).sum()[0]/float(df.shape[0])) + '\n')

        fout.close()


parser = argparse.ArgumentParser(description='Compute and plot BaseQuality statistics from FastQ or BAM file')
parser.add_argument(
    '--fastq',
    type=str,
    required=False,
    nargs='+',
    default=None,
    help='Space-delimited list of on or more Fastq files'
)
parser.add_argument(
    '--gz',
    action='store_true',
    required=False,
    help='If fastq files are compressed (.gz)'
)
parser.add_argument(
    '--bam',
    type=str,
    required=False,
    default=None,
    help='BAM file'
)
parser.add_argument(
    '--sample',
    type=str,
    required=False,
    default='',
    help='Sample name'
)
parser.add_argument(
    '--prefix',
    type=str,
    required=True,
    help='Output prefix name'
)
args = parser.parse_args()

if args.bam is not None and args.fastq is not None:
    print 'Supply only BAM or FastQ file!'
    sys.exit()

main(args=args)

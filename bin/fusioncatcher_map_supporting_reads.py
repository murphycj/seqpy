import sys
import pandas
import numpy as np
import os
import argparse

def align_with_star(star,fastq_file,genome_dir,bam_prefix):

    cmd = star + ' ' + \
    '--outSAMtype BAM SortedByCoordinate ' + \
    '--runThreadN 1 ' + \
    '--genomeDir ' + genome_dir + ' ' + \
    '--readFilesIn ' + fastq_file + ' ' + \
    '--outFileNamePrefix ' + bam_prefix
    os.system(cmd)

    os.system(args.samtools + ' index ' + bam_prefix + 'Aligned.sortedByCoord.out.bam')

def bam_to_fastq(bam_file):
    os.system(args.samtools + ' bam2fq tmp.bam > tmp.fastq')
    os.system('gzip tmp.fastq')
    return 'tmp.fastq.gz'

def write_fasta(fusion,seq):
    fasta_name = fusion + '/' + fusion + '.fa'
    fout = open(fasta_name,'w')
    fout.write('>' + fusion + '\n')
    fout.write(seq + '\n')
    fout.close()

    return fasta_name

def make_star_reference(star,fasta,genome_dir):
    os.system(
        star + ' ' +
        '--runThreadN 1 ' + \
        '--runMode genomeGenerate ' + \
        '--genomeDir ' + genome_dir + ' ' + \
        '--genomeFastaFiles ' + fasta
    )

def main(args):

    data = dict()
    input_data = dict()
    read_counts = dict()

    genome_dir = args.fusion + '/star'
    bam_prefix = args.fusion + '/' + args.fusion + '.'

    if not os.path.exists(args.fusion):
        os.mkdir(args.fusion)
    if not os.path.exists(genome_dir):
        os.mkdir(genome_dir)

    #fasta_name = write_fasta(args.fusion,args.seq)
    make_star_reference(args.star,args.fasta,genome_dir)

    align_with_star(star=args.star,fastq_file=args.fq_file,genome_dir=genome_dir,bam_prefix=bam_prefix)

parser = argparse.ArgumentParser(description='Counts fusion supporting and non-supporting reads.')
parser.add_argument(
    '--fusion',
    type=str,
    required=True,
    help='prefix of anaylsis'
)
parser.add_argument(
    '--fasta',
    type=str,
    required=True,
    help='Fused seqeuence'
)
parser.add_argument(
    '--fq_file',
    type=str,
    required=True,
    help='Fastq file containing supporting reads'
)
parser.add_argument(
    '--star',
    type=str,
    required=True,
    help='Path to STAR'
)
parser.add_argument(
    '--samtools',
    type=str,
    required=True,
    help='Path to samtools'
)
args = parser.parse_args()

main(args=args)

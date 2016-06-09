"""
This script will take the STAR chimeric alignment and count the reads that
do support and do not support the give gene fusion

Does not support PE reads for now
"""

import sys
import pandas
import numpy as np
import os
import argparse
import pysam

def get_and_print_reads(bam,bamchimeric,args):
    bam_file = pysam.Samfile(bam,'rb')

    if os.path.splitext(bamchimeric)[1]=='.sam':

        print 'Chimeric reads are in a SAM file, converting to BAM...'

        os.system(args.samtools + ' view -bS ' + bamchimeric + ' > chimeric.bam')
        os.system(args.samtools + ' sort chimeric.bam chimeric.sorted')
        os.system(args.samtools + ' index chimeric.sorted.bam')

        bam_chimeric_file = pysam.Samfile('chimeric.sorted.bam','rb')

    else:
        bam_chimeric_file = pysam.Samfile(bamchimeric,'rb')

    bam_out = pysam.Samfile('tmp.bam','wb',header=bam_file.header)

    reads = []

    #minus one because pysam uses 0-based

    for read in bam_file.fetch(args.gene1chrom,int(args.gene1start)-1,int(args.gene1end)-1):
        reads.append(read.qname)
        bam_out.write(read)
    for read in bam_file.fetch(args.gene2chrom,int(args.gene2start)-1,int(args.gene2end)-1):
        reads.append(read.qname)
        bam_out.write(read)

    for read in bam_chimeric_file.fetch(args.gene1chrom,int(args.gene1start)-1,int(args.gene1end)-1):
        if read.qname not in reads:
            bam_out.write(read)
    for read in bam_chimeric_file.fetch(args.gene2chrom,int(args.gene2start)-1,int(args.gene2end)-1):
        if read.qname not in reads:
            bam_out.write(read)

    bam_file.close()
    bam_chimeric_file.close()
    bam_out.close()

    os.system(args.samtools + ' index tmp.bam')

    return 'tmp.bam'

def align_with_star(fastq_file,sample_out_dir):
    out_bam = './' + sample_out_dir + '/' + sample_out_dir + '.Aligned.sortedByCoord.out.bam'
    os.system(
        args.STAR + ' ' + \
        '--outSAMtype BAM SortedByCoordinate ' + \
        '--runThreadN 1 ' + \
        '--genomeDir ' + args.STARref + ' ' + \
        '--readFilesCommand zcat ' + \
        '--readFilesIn ' + fastq_file + ' ' + \
        '--outFileNamePrefix ./' + sample_out_dir + '/' + sample_out_dir + '.'
    )

    os.system(args.samtools + ' index ' + out_bam)

    return out_bam

def count_alignments(new_bam,args):
    #some reads are said to align to the fusion sequence, but the portion
    # of the sequence that maps after/before the junction are soft clipped
    #will count these sequences as mapping to the gene with the most
    #matches

    fusion_name = args.gene1name + '-' + args.gene2name

    counts = {
        args.gene1name:0,
        args.gene2name:0,
        fusion_name:0
    }

    bam_file = pysam.Samfile(new_bam,'rb')

    r_names = map(lambda x: x['SN'], bam_file.header['SQ'])

    for alignment in bam_file.fetch():

        rname = r_names[alignment.rname]

        if rname==args.gene1name or rname==args.gene2name:
            counts[rname]+=1
        else:

            #is it mostly in gene1 or gene2?
            gene1 = False
            gene2 = False

            if (args.junctionPosition - alignment.pos) > (alignment.qlen/2):
                #is the portion in gene2 soft-clipped?

                length_not_soft_clipped = 0
                for c in alignment.cigar:
                    if c[0]==0:
                        length_not_soft_clipped += c[1]

                stop_of_none_clipped = alignment.pos + length_not_soft_clipped

                if stop_of_none_clipped!=args.junctionPosition:
                    counts[fusion_name]+=1
                else:
                    counts[args.gene1name]+=1
            else:
                #is the portion in gene1 soft-clipped?

                counts[fusion_name]+=1

    bam_file.close()
    return counts

def write_results(read_counts,args):
    fout = open(args.out,'w')

    fusion_name = args.gene1name + '-' + args.gene2name

    fout.write('sample,number fusion supporting reads,number non-fusion supporting reads (' + args.gene1name + '),number non-fusion supporting reads (' + args.gene2name + ')\n')
    for i,j in read_counts.items():
        fout.write(i + ',')
        fout.write(str(j[fusion_name]) + ',')
        fout.write(str(j[args.gene1name]) + ',')
        fout.write(str(j[args.gene2name]) + '\n')
    fout.close()

def bam_to_fastq(bam_file):
    os.system(args.samtools + ' bam2fq tmp.bam > tmp.fastq')
    os.system('gzip tmp.fastq')
    return 'tmp.fastq.gz'

def main(args):

    data = dict()
    input_data = dict()
    read_counts = dict()
    samples = args.samples
    bams = args.bam
    bamChimerics = args.bamChimeric
    for i in range(0,len(samples)):
        data[samples[i]] = [
            bams[i],
            bamChimerics[i]
        ]

        fusion_name = args.gene1name + '-' + args.gene2name

        read_counts[samples[i]] = {
            args.gene1name:0,
            args.gene1name:0,
            fusion_name: 0
        }


    for s in data.keys():

        sample_out_dir = s + '.' + args.gene1name + '-' + args.gene2name

        if not os.path.exists(sample_out_dir):
            os.mkdir(sample_out_dir)

        fusion_name = args.gene1name + '-' + args.gene2name

        #get the reads

        bam_file = get_and_print_reads(bam=data[s][0],bamchimeric=data[s][1],args=args)

        #convert to fastq

        fastq_file = bam_to_fastq(bam_file=bam_file)

        #realign

        new_bam = align_with_star(fastq_file=fastq_file,sample_out_dir=sample_out_dir)

        #count alignments

        counts = count_alignments(new_bam=new_bam,args=args)

        read_counts[s][args.gene1name] = counts[args.gene1name]
        read_counts[s][args.gene2name] = counts[args.gene2name]
        read_counts[s][fusion_name] = counts[fusion_name]

        #clean up
        os.system('rm tmp.*')

        if os.path.exists('chimeric.bam'):
            os.system('rm chimeric.bam')
            os.system('rm chimeric.sorted.bam')
            os.system('rm chimeric.sorted.bam.bai')

    write_results(read_counts=read_counts,args=args)


parser = argparse.ArgumentParser(description='Counts fusion supporting and non-supporting reads.')
parser.add_argument(
    '--samples',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of sample names.'
)
parser.add_argument(
    '--bam',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of BAM files'
)
parser.add_argument(
    '--bamChimeric',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of chimeric alignments from STAR in BAM files'
)
parser.add_argument(
    '--gene1name',
    type=str,
    required=True,
    help='Name of gene 1'
)
parser.add_argument(
    '--gene1chrom',
    type=str,
    required=True,
    help='Chromosome of gene 1'
)
parser.add_argument(
    '--gene1start',
    type=str,
    required=True,
    help='Start position of gene 1 junction (1-based)'
)
parser.add_argument(
    '--gene1end',
    type=str,
    required=True,
    help='End position of gene 1 junction (1-based)'
)
parser.add_argument(
    '--gene2name',
    type=str,
    required=True,
    help='Name of gene 2'
)
parser.add_argument(
    '--gene2chrom',
    type=str,
    required=True,
    help='Chromosome of gene 2'
)
parser.add_argument(
    '--gene2start',
    type=str,
    required=True,
    help='Start position of gene 2 junction (1-based)'
)
parser.add_argument(
    '--gene2end',
    type=str,
    required=True,
    help='End position of gene 2 junction (1-based)'
)
parser.add_argument(
    '--junctionPosition',
    type=int,
    required=True,
    help='Position of the junction in the fusion sequence (1-based)'
)
parser.add_argument(
    '--STAR',
    type=str,
    required=True,
    help='Path to STAR'
)
parser.add_argument(
    '--STARref',
    type=str,
    required=True,
    help='Reference for STAR'
)
parser.add_argument(
    '--samtools',
    type=str,
    required=True,
    help='Path to samtools'
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file to save results'
)
args = parser.parse_args()

main(args=args)

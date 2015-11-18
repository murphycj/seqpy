import multiprocessing
import sys
import pandas
import glob
import os
import argparse

def job(filename):
    os.system('./' + filename)

def batch_run(args):

    dirs = glob.glob(args.meta + '/' + args.dirpattern + '*')
    pool = multiprocessing.Pool(args.p)

    for f in dirs:
        if not os.path.isdir(f):
            continue

        sample = os.path.split(f)[1]
        bam = glob.glob(f + '/*' + args.pattern)
        if len(bam)!=1:
            print 'Problem finding bam file'
            continue

        if not os.path.exists(f + '/HTSeqCount/'):
            os.mkdir(f + '/HTSeqCount')
        outfile = f + '/HTSeqCount/' + sample + '.' + args.pattern + '.count'

        filename = sample + '/HTSeqCount/' + sample + '.htseq.count.txt'
        fout = open(filename,'w')

        fout.write(
            'python ~/.local/bin/htseq-count ' + \
            '-s no ' + \
            '-f bam ' + \
             bam[0] + ' ' + \
             args.gff + ' ' + \
             '> ' + outfile + '\n'
        )

        fout.close()
        os.system('chmod 755 ' + filename)

        pool.apply_async(job,args=(filename,))
    pool.close()
    pool.join()


parser = argparse.ArgumentParser(description='Aggregates the FPKM values from cufflinks output')
parser.add_argument('--meta',type=str,help='Meta-directory where each sub-directory contains bam for each sample',required=True)
parser.add_argument('--gff',type=str,help='gff file',required=True)
parser.add_argument('--pattern',type=str,help='BAM file patter',required=True)
parser.add_argument('--dirpattern',type=str,help='BAM file patter',required=True)
parser.add_argument('--p',type=int,help='Number of threads',default=6)
args = parser.parse_args()

args.meta = os.path.abspath(args.meta)
print args.meta
batch_run(args=args)

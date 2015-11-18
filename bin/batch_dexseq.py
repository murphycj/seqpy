import multiprocessing
import sys
import pandas
import glob
import os
import argparse

def job(filename):
    os.system('./' + filename)

def batch_run(args):

    dirs = glob.glob(args.meta + '/' + 'HL*')
    pool = multiprocessing.Pool(args.p)

    for f in dirs:
        if not os.path.isdir(f):
            continue

        sample = os.path.split(f)[1]
        print sample
        continue
        bam = glob.glob(f + '/' + args.pattern)
        if len(bam)!=1:
            print 'BAM file not found, ' + sample
            continue
        bam = bam[0]
        count_file = f + '/' + sample + '.Aligned.sortedByCoord.out.count'
        filename = sample + '/' + sample + '_dexseq_count.sh'

        fout = open(filename,'w')
        fout.write(
            '/home/chm2059/chm2059/lib/python2.7.8/bin/python ' + \
            '/home/chm2059/chm2059/lib/DEXSeq/python_scripts/dexseq_count.py ' + \
            '-f bam ' + \
            '-s no ' + \
            '-r pos ' + \
            args.gff + \
            ' ' + \
            bam + \
            ' ' + \
            count_file + '\n'
        )
        fout.close()
        os.system('chmod 755 ' + filename)

        pool.apply_async(job,args=(filename,))
    #pool.close()
    #pool.join()


parser = argparse.ArgumentParser(description='Parallel running of DEXseq_count.py')
parser.add_argument('--meta',type=str,help='Meta-directory where each sub-directory contains bam for each sample',required=True)
parser.add_argument('--pattern',type=str,help="String pattern used to search for BAM files",required=True)
parser.add_argument('--gff',type=str,help='gff file',required=True)
parser.add_argument('--p',type=int,help='Number of threads',default=6)
args = parser.parse_args()

args.meta = os.path.abspath(args.meta)
print args.meta
batch_run(args=args)

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

        if not os.path.exists(f + '/DEXseq/'):
            os.mkdir(f + '/DEXseq')
        outfile = f + '/DEXseq/' + sample + '.' + args.pattern + '.count'

        filename = sample + '/DEXseq/' + sample + '_dexseq_count.sh'
        tmp_bam = sample + '/DEXseq/tmp.bam'

        if args.addNHfield:
            fout = open(filename,'w')
            fout.write(
                '/home/chm2059/chm2059/lib/python2.7.8/bin/python ' + \
                '~/chm2059/lib/seqpy/bin/prepare_bam_for_dexseq.py ' + \
                '--bam ' + bam[0] + ' ' + \
                '--out ' + tmp_bam + ' ' + \
                '--XTtag\n'
            )
            fout.write(
                '/home/chm2059/chm2059/lib/python2.7.8/bin/python ' + \
                '/home/chm2059/chm2059/lib/DEXSeq/python_scripts/dexseq_count.py ' + \
                '-f bam ' + \
                '-s no ' + \
                '-r pos ' + \
                args.gff + \
                ' ' + \
                tmp_bam + \
                ' ' + \
                outfile + '\n'
            )
            fout.close()
            os.system('chmod 755 ' + filename)


        else:
            fout = open(filename,'w')
            fout.write(
                '/home/chm2059/chm2059/lib/python2.7.8/bin/python ' + \
                '/home/chm2059/chm2059/lib/DEXSeq/python_scripts/dexseq_count.py ' + \
                '-f bam ' + \
                '-s no ' + \
                '-r pos ' + \
                args.gff + \
                ' ' + \
                bam[0] + \
                ' ' + \
                outfile + '\n'
            )
            fout.close()
            os.system('chmod 755 ' + filename)

        pool.apply_async(job,args=(filename,))
    pool.close()
    pool.join()


parser = argparse.ArgumentParser(description='Parallel running of DEXseq_count.py')
parser.add_argument('--meta',type=str,help='Meta-directory where each sub-directory contains bam for each sample',required=True)
parser.add_argument('--pattern',type=str,help='BAM file pattern',required=True)
parser.add_argument('--dirpattern',type=str,help='directory file pattern',required=True)
parser.add_argument('--gff',type=str,help='gff file',required=True)
parser.add_argument('--p',type=int,help='Number of threads',default=6)
parser.add_argument('--addNHfield',action='store_true',help='Whether the BAM file needs the NH:i:1 flag added to the reads.')
args = parser.parse_args()

args.meta = os.path.abspath(args.meta)
print args.meta

batch_run(args=args)

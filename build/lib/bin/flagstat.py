import multiprocessing
import argparse
import pandas
from crimson import flagstat
import glob
import os

def job(filename):
    os.system(filename)

def run_flagstat(args):
    dirs = glob.glob(args.dir + '/*')
    pool = multiprocessing.Pool(args.p)

    for f in dirs:
        if not os.path.isdir(f):
            continue

        sample = os.path.split(f)[1]

        bam = glob.glob(f + '/*' + args.bam)
        if len(bam)!=1:
            print 'Found the following for ' + sample
            continue

        filename = f + '/' + sample + '.' + args.bam + '.flagstat.sh'
        fout = open(filename,'w')
        fout.write(
            '/home/chm2059/chm2059/lib/samtools-1.2/samtools flagstat ' + \
            bam[0] + \
             '> ' + bam[0] + '.flagstat\n'
        )
        fout.close()
        os.system('chmod 755 ' + filename)
        pool.apply_async(job,args=(filename,))
    pool.close()
    pool.join()


def summarize_flagstat(args):
    dirs = glob.glob(args.dir + '/*')
    data = {}

    for f in dirs:
        if not os.path.isdir(f):
            continue
        sample = os.path.split(f)[1]
        files = glob.glob(f + '/*' + args.bam + '.flagstat')
        if len(files)==0:
            continue
        data[sample] = flagstat.parse(files[0])['pass_qc']
    results = pandas.DataFrame(data).T
    results.to_csv(args.out)

parser = argparse.ArgumentParser(description='Run and or summarize samtools flagstat results')
parser.add_argument('--dir',type=str,help='Meta-directory where each sub-directory contains the BAM file',required=True)
parser.add_argument('--bam',type=str,help='BAM file name pattern (e.g. --bam sorted.bam)',required=True)
parser.add_argument('--summarize',action='store_true',help='Summarize flagstat',default=False)
parser.add_argument('--run',action='store_true',help='Summarize flagstat',default=False)
parser.add_argument('--out',type=str,help='Output file name',required=True)
parser.add_argument('--p',type=int,help='Number processors',required=True)
args = parser.parse_args()

if args.run:
    run_flagstat(args=args)

if args.summarize:
    summarize_flagstat(args=args)

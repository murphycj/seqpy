import numpy as np
import argparse
import os
import pandas
import mouse2human


CHIP = 'gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip'
MSIGDB_ROOT = '/Users/charlesmurphy/Desktop/Research/data/msigdb/'

def create_rnk(args,data):

    gct_file = os.path.abspath(args.outDir) + '/' + args.outDir + '.rnk'

    if not args.human:
        c = mouse2human.Conversion()

        data = c.convert_mouse2human(data=data,remove_duplicates=True)

    result_data = pandas.DataFrame(columns=['name','rank'],index=data.index)


    if args.limmavoom:
        if args.ranking == 'statistic':
            result_data['name'] = data.index
            result_data['rank'] = data['t'].tolist()
            result_data = result_data[~pandas.isnull(result_data['name'])]
            result_data = result_data[~pandas.isnull(result_data['rank'])]
        elif args.ranking == 'logpvalFC':
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['P.Value'])*data['logFC']).tolist()
            result_data = result_data[~pandas.isnull(result_data['name'])]
            result_data = result_data[~pandas.isnull(result_data['rank'])]
    else:
        if args.ranking == 'statistic':
            result_data['name'] = data.index
            result_data['rank'] = data['stat'].tolist()
            result_data = result_data[~pandas.isnull(result_data['name'])]
            result_data = result_data[~pandas.isnull(result_data['rank'])]
        elif args.ranking == 'logpvalFC':
            result_data['name'] = data.index
            result_data['rank'] = (-np.log2(data['pvalue'])*data['log2FoldChange']).tolist()
            result_data = result_data[~pandas.isnull(result_data['name'])]
            result_data = result_data[~pandas.isnull(result_data['rank'])]

    rnk_file = args.outDir + '/' + args.outDir + '_' + args.msigdb.replace('.','_') + '.rnk'
    result_data.to_csv(rnk_file,index=False,header=False, sep='\t')

    return rnk_file

def write_script(args,data,rnk_file):


    #write script
    if args.cluster:
        filename = os.path.abspath(args.outDir) + '/' + args.outDir + '-' + args.msigdb + '.sh'
        fout = open(filename,'w')
        fout.write('#! /bin/bash -l\n')
        fout.write('#$ -j y\n')
        fout.write('#$ -o /home/chm2059/igv-batch.log\n')
        fout.write('#$ -l h_rt=1:00:00\n')
        fout.write('#$ -pe smp 1\n')
        fout.write('#$ -l h_vmem=8G\n')
        fout.write('#$ -l os=rhel5.4|rhel6.3\n')
        fout.write('#$ -R y\n')
        fout.write('# fpkm : ' + args.fpkm + '\n')
        fout.write('# GSEA : '  + os.path.abspath(args.gsea) + '\n')
        fout.write('# group1 : ' + args.group1 + '\n')
        fout.write('# group2 : ' + args.group2 + '\n')
        fout.write('# msigdb : ' + args.msigdb + '\n')
        fout.write('# fpkmThreshold : ' + str(args.fpkmThreshold) + '\n')
        fout.write('# sampleNonZero : ' + str(args.sampleNonZero) + '\n')
        fout.write('# outDir : ' + args.outDir + '\n\n')
        fout.write('rsync ' + gct_file + ' $TMPDIR\n')
        fout.write('rsync ' + cls_file + ' $TMPDIR\n')
        fout.write('cd $TMPDIR\n')
        fout.write(
            'java -cp' + os.path.abspath(args.gsea) + '-Xmx2048m xtools.gsea.GseaPreranked ' + \
            '-res ' + os.path.split(gct_file)[1] + ' ' + \
            '-cls ' + os.path.split(cls_file)[1] + '#' + phenotypes[0] + '_versus_' + phenotypes[1] + \
            '-gmx ' + MSIGDB_ROOT + args.msigdb + \
            ' -collapse false ' + \
            '-mode Max_probe ' + \
            '-norm meandiv ' + \
            '-nperm 1000 ' + \
            '-permute phenotype ' + \
            '-rnd_type no_balance ' + \
            '-scoring_scheme weighted ' + \
            '-rpt_label ' + args.outDir + '-' + args.msigdb + ' '\
            '-metric Signal2Noise ' + \
            '-sort real ' + \
            '-order descending ' + \
            '-include_only_symbols true ' + \
            '-make_sets true ' + \
            '-median false ' + \
            '-num 100 ' + \
            '-plot_top_x 40 ' + \
            '-rnd_seed timestamp ' + \
            '-save_rnd_lists false ' + \
            '-set_max 500 ' + \
            '-set_min 15 ' + \
            '-zip_report true ' + \
            '-out $TMPDIR/' + \
            ' -gui false\n'
        )
        fout.write('rm *gct\n')
        fout.write('rm *cls\n')
        fout.write('rsync ' + args.outDir + '* ' + os.path.abspath(os.outDir) + '\n')
        fout.close()
        os.system('chmod 755 ' + filename)
    else:
        filename = os.path.abspath(args.outDir) + '/' + args.outDir + '-' + args.msigdb + '.sh'
        fout = open(filename,'w')
        fout.write('#! /bin/bash -l\n')
        fout.write('# GSEA : '  + os.path.abspath(args.gsea) + '\n')
        fout.write('# msigdb : ' + args.msigdb + '\n')
        fout.write('# outDir : ' + args.outDir + '\n\n')
        fout.write(
            'java -cp ' + os.path.abspath(args.gsea) + ' -Xmx2048m xtools.gsea.GseaPreranked \\\n' + \
            ' -rnk ' + rnk_file + ' \\\n' + \
            ' -gmx ' + MSIGDB_ROOT + args.msigdb + ' \\\n' + \
            ' -collapse false \\\n' + \
            ' -mode Max_probe \\\n' + \
            ' -norm meandiv \\\n' + \
            ' -nperm 1000 \\\n' + \
            ' -scoring_scheme weighted \\\n' + \
            ' -rpt_label ' + args.outDir + '_' + args.msigdb + ' \\\n'\
            ' -include_only_symbols true \\\n' + \
            ' -make_sets true \\\n' + \
            ' -plot_top_x 40 \\\n' + \
            ' -rnd_seed timestamp \\\n' + \
            ' -set_max 500 \\\n' + \
            ' -set_min 15 \\\n' + \
            ' -zip_report true \\\n' + \
            ' -out ' + os.path.abspath(args.outDir) + ' \\\n' +\
            ' -gui false\n'
        )
        fout.close()
        os.system('chmod 755 ' + filename)

def main(args):

    if args.cluster:
        MSIGDB_ROOT = '/home/chm2059/chm2059/data/msigdb/'

    data = pandas.read_csv(args.deseq2, header=0, index_col=0)
    data = data.dropna()

    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)

    rnk_file = create_rnk(args=args,data=data)

    write_script(args=args,data=data,rnk_file=rnk_file)


parser = argparse.ArgumentParser(description='Creates GSEA scripts to run')
parser.add_argument('--deseq2',type=str,help='csv file containing deseq2 results',required=True)
parser.add_argument('--gsea',type=str,help='Path to GSEA JAR',required=True)
parser.add_argument('--msigdb',type=str,help='The MSigDB to use (e.g. c2.cp.v5.0.symbols.gmt)',required=True)
parser.add_argument('--cluster',action='store_true',help='Write scripts for running as jobs on cluster (default False)',required=False, default=False)
parser.add_argument('--ranking',type=str,help='Ranking statistic to use (statistic, logpvalFC)',default='statistic')
parser.add_argument('--limmavoom',action='store_true',help='If results are from limmavoom')
parser.add_argument('--human',action='store_true',help='If the gene symbols are human',default=False)
parser.add_argument('--outDir',type=str,help='Output directory',required=True)
args = parser.parse_args()

main(args=args)

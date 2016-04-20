import argparse
import os
import pandas

CHIP = 'gseaftp.broadinstitute.org://pub/gsea/annotations/GENE_SYMBOL.chip'
MSIGDB_ROOT = '/Users/charlesmurphy/Desktop/Research/data/msigdb/'

def create_gct(args,fpkm):

    gct_file = os.path.abspath(args.outDir) + '/' + args.outDir + '.gct'

    samples = fpkm.columns.tolist()
    genes = fpkm.index.tolist()
    descriptions = ['NA']*fpkm.shape[0]
    temp = ['NAME','DESCRIPTION'] + samples

    gct = pandas.DataFrame(columns=temp, index=fpkm.index)
    gct['NAME'] = genes
    gct['DESCRIPTION'] = descriptions
    for i in samples:
      gct[i] = fpkm[i]
    gct.to_csv(os.path.abspath(args.outDir) + '/temp.txt',sep='\t',index=False)

    gct = open(os.path.abspath(args.outDir) + '/temp.txt', 'r').read()
    fout = open(gct_file, 'w')
    fout.write('#1.2\n')
    fout.write(str(fpkm.shape[0]) + '\t' + str(len(samples)) + '\n')
    fout.write(gct)
    fout.close()
    os.system('rm ' + os.path.abspath(args.outDir) + '/temp.txt')

    return gct_file

def create_cls(args, fpkm, group1, group2,phenotypes):

    cls_file = os.path.abspath(args.outDir) + '/' + args.outDir + '.cls'

    classes = [phenotypes[0],phenotypes[1]]
    class_labels = [phenotypes[0]]*len(group1) + [phenotypes[1]]*len(group2)

    fout = open(cls_file, 'w')
    fout.write(str(len(class_labels)) + '\t' + str(len(classes)) + '\t1\n')
    fout.write('#\t' + '\t'.join(classes) + '\n' + class_labels[0])
    for i in range(1, len(class_labels)):
      fout.write('\t' + class_labels[i])
    fout.write('\n')
    fout.close()

    return cls_file

def write_script(args,fpkm,gct_file,cls_file,phenotypes):


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
            'java -cp' + os.path.abspath(args.gsea) + '-Xmx2048m xtools.gsea.Gsea ' + \
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
            '-plot_top_x 20 ' + \
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
        fout.write('# fpkm : ' + args.fpkm + '\n')
        fout.write('# GSEA : '  + os.path.abspath(args.gsea) + '\n')
        fout.write('# group1 : ' + args.group1 + '\n')
        fout.write('# group2 : ' + args.group2 + '\n')
        fout.write('# msigdb : ' + args.msigdb + '\n')
        fout.write('# fpkmThreshold : ' + str(args.fpkmThreshold) + '\n')
        fout.write('# sampleNonZero : ' + str(args.sampleNonZero) + '\n')
        fout.write('# outDir : ' + args.outDir + '\n\n')
        fout.write(
            'java -cp ' + os.path.abspath(args.gsea) + ' -Xmx2048m xtools.gsea.Gsea \\\n' + \
            ' -res ' + gct_file + ' \\\n' + \
            ' -cls ' + cls_file + '#' + phenotypes[0] + '_versus_' + phenotypes[1] + ' \\\n' + \
            ' -gmx ' + MSIGDB_ROOT + args.msigdb + ' \\\n' + \
            ' -collapse false \\\n' + \
            ' -mode Max_probe \\\n' + \
            ' -norm meandiv \\\n' + \
            ' -nperm 1000 \\\n' + \
            ' -permute phenotype \\\n' + \
            ' -rnd_type no_balance \\\n' + \
            ' -scoring_scheme weighted \\\n' + \
            ' -rpt_label ' + args.outDir + '_' + args.msigdb + ' \\\n'\
            ' -metric Signal2Noise \\\n' + \
            ' -sort real \\\n' + \
            ' -order descending \\\n' + \
            ' -include_only_symbols true \\\n' + \
            ' -make_sets true \\\n' + \
            ' -median false \\\n' + \
            ' -num 100 \\\n' + \
            ' -plot_top_x 20 \\\n' + \
            ' -rnd_seed timestamp \\\n' + \
            ' -save_rnd_lists false \\\n' + \
            ' -set_max 500 \\\n' + \
            ' -set_min 15 \\\n' + \
            ' -zip_report true \\\n' + \
            ' -out ' + os.path.abspath(args.outDir) + ' \\\n' +\
            ' -gui false\n'
        )
        fout.close()
        os.system('chmod 755 ' + filename)

def main(args):

    group1 = args.group1.split(',')
    group2 = args.group2.split(',')

    phenotypes = args.phenotypes.split(',')
    assert len(phenotypes)==2, 'supply only two phenotypes'

    if args.cluster:
        MSIGDB_ROOT = '/home/chm2059/chm2059/data/msigdb/'

    fpkm = pandas.read_csv(args.fpkm,index_col=0)

    for i in group1 + group2:
        assert i in fpkm.columns, i + " not in fpkm"

    fpkm = fpkm[group1 + group2]
    fpkm[fpkm<=args.fpkmThreshold]=0
    fpkm = fpkm[fpkm.sum(axis=1)>=args.sampleNonZero]

    if not os.path.exists(args.outDir):
        os.mkdir(args.outDir)

    gct_file = create_gct(args=args,fpkm=fpkm)

    cls_file = create_cls(args=args,fpkm=fpkm, group1=group1, group2=group2, phenotypes=phenotypes)

    write_script(args=args,fpkm=fpkm,gct_file=gct_file,cls_file=cls_file, phenotypes=phenotypes)


parser = argparse.ArgumentParser(description='Creates GSEA scripts to run')
parser.add_argument('--fpkm',type=str,help='csv file containing FPKM values (using human gene symbols)',required=True)
parser.add_argument('--gsea',type=str,help='Path to GSEA JAR',required=True)
parser.add_argument('--group1',type=str,help='Comma separated group 1 samples',required=True)
parser.add_argument('--group2',type=str,help='Comma separated group 2 samples',required=True)
parser.add_argument('--phenotypes',type=str,help='Comma separated name for the groups/phenotype (e.g. WT,MUT), same order and group1 and group2',required=True)
parser.add_argument('--msigdb',type=str,help='The MSigDB to use (e.g. c2.cp.v5.0.symbols.gmt)',required=True)
parser.add_argument('--fpkmThreshold',type=float,help='FPKM threshold',required=False,default=0.1)
parser.add_argument('--sampleNonZero',type=int,help='Minimum number of samples with non-zeros expression (defualt 1)',required=False,default=1)
parser.add_argument('--cluster',action='store_true',help='Write scripts for running as jobs on cluster (default False)',required=False, default=False)
parser.add_argument('--outDir',type=str,help='Output directory',required=True)
args = parser.parse_args()

main(args=args)

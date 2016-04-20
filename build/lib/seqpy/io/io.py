"""
rnaseq-lib.py
test
Provides functions that I used to process and plot the data

"""

import os
import __main__
import logging
import git
import datetime
import atexit
import glob
import shutil
import imp
import sys
import random

import pandas

try:
  import matplotlib.pyplot as plt
  import matplotlib
  from matplotlib.mlab import PCA
  from pylab import get_cmap
except:
  pass

try:
  import scipy
  from scipy.spatial.distance import pdist, squareform
  from scipy.cluster.hierarchy import linkage, dendrogram
except:
  pass

import numpy as np
import xml
import xml.etree.ElementTree as ET

LIBPATH = os.path.dirname(os.path.realpath(__file__))
comparisons = imp.load_source('comparisons',LIBPATH + '/comparisons.py')
unittests = imp.load_source('unittests',LIBPATH + '/unittests.py')
mutations = imp.load_source('mutations',LIBPATH + '/mutations.py')

log = logging

SAMPLE_INFO = pandas.read_csv(LIBPATH + '/samples/samples.csv')


CLUSTER_SRCDIR = '/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/ole2001/'
CLUSTER_WORKDIR = '/zenodotus/dat02/elemento_lab_scratch/oelab_scratch_scratch007/chm2059/elementoCantley_2014_9_29/results/9_30_2014-copydata/'

DIRS = ['HuiLiu_2014_05_12',
        'HuiLiu_2014_06_24',
        'HuiLiu_2014_07_24',
        'HuiLiu_2014_08_06',
        'HuiLiu_2014_09_05']


def _get_time():
  today = str(datetime.datetime.today())
  today = today.replace(' ', '-')
  today = today.replace(':', '.')
  return today

def _compute_mouse_human_homolog(mouse_list, idType='EntrezGene ID', conversion_file=''):
  """
  Convert a list of mouse gene symbols to human symbols using the
  MGI homolog database:
  ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt

  mouse_list: list of mouse identifiers
  idType: type of ID to use. Available: Symbol, EntrezGene ID, OMIM Gene ID, HGNC ID

  """

  MOUSE_HUMAN_HOMOLOGS = pandas.read_table('../../data/10_21_14-HOM_MouseHumanSequence.rpt',sep='\t')
  mouse_homologs = MOUSE_HUMAN_HOMOLOGS[MOUSE_HUMAN_HOMOLOGS['Common Organism Name']=='mouse, laboratory']
  human_homologs = MOUSE_HUMAN_HOMOLOGS[MOUSE_HUMAN_HOMOLOGS['Common Organism Name']=='human']

  homologs = []
  filtered_mouse_list = []

  count = 0
  for gene in mouse_list:
    temp = mouse_homologs[mouse_homologs['Symbol']==gene]

    #check to see if there is an entry for the mouse symbol
    if temp.shape[0]==0:
      homologs.append('NA')
      filtered_mouse_list.append('NA')
      count += 1
      continue

    temp = temp['HomoloGene ID'].tolist()[0]
    temp = human_homologs[human_homologs['HomoloGene ID']==temp]

    #check to see if there is a human homolog
    if temp.shape[0]!=0:
      homologs.append(temp[idType].tolist()[0])
      filtered_mouse_list.append(gene)
    else:
      count += 1
      homologs.append('NA')
      filtered_mouse_list.append('NA')
  log.info('WARNING: could not find either mouse entry or human homolog for ' + str(count) + ' mouse ids.')

  save_data = pandas.DataFrame(index=range(0,len(homologs)),columns=['human_homolog','mouse_ids'])
  save_data['human_homolog'] = homologs
  save_data['mouse_ids'] = filtered_mouse_list
  save_data.to_csv(conversion_file)

def _compute_human_to_mouse_homolog(human_list, idType='EntrezGene ID', conversion_file=''):
  """
  Convert a list of human gene symbols to mouse symbols using the
  MGI homolog database:
  ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt

  mouse_list: list of mouse identifiers
  idType: type of ID to use. Available: Symbol, EntrezGene ID, OMIM Gene ID, HGNC ID

  """

  MOUSE_HUMAN_HOMOLOGS = pandas.read_table('../../data/10_21_14-HOM_MouseHumanSequence.rpt',sep='\t')
  mouse_homologs = MOUSE_HUMAN_HOMOLOGS[MOUSE_HUMAN_HOMOLOGS['Common Organism Name']=='mouse, laboratory']
  human_homologs = MOUSE_HUMAN_HOMOLOGS[MOUSE_HUMAN_HOMOLOGS['Common Organism Name']=='human']

  homologs = []
  filtered_human_list = []

  count = 0
  for gene in human_list:

    temp = human_homologs[human_homologs[idType]==gene]

    #check to see if there is an entry for the mouse symbol
    if temp.shape[0]==0:
      homologs.append('NA')
      filtered_human_list.append('NA')
      count += 1
      continue

    temp = temp['HomoloGene ID'].tolist()[0]
    temp = mouse_homologs[mouse_homologs['HomoloGene ID']==temp]

    #check to see if there is a human homolog
    if temp.shape[0]!=0:
      homologs.append(temp[idType].tolist()[0])
      filtered_human_list.append(gene)
    else:
      count += 1
      homologs.append('NA')
      filtered_human_list.append('NA')
  log.info('WARNING: could not find either mouse entry or human homolog for ' + str(count) + ' human ids.')

  save_data = pandas.DataFrame(index=range(0,len(homologs)),columns=['mouse_homolog','human_ids'])
  save_data['mouse_homolog'] = homologs
  save_data['human_ids'] = filtered_human_list
  save_data.to_csv(conversion_file)


def setup_log(log):
  """
  Sets up the log file everytime the library is loaded, save the log
  file to the ./log directory
  """

  #get current git version
  repo = git.Repo(LIBPATH)
  version = str(repo.commit('master').hexsha)
  today = _get_time()
  #create log file
  if not os.path.exists('log'):
      os.mkdir('log')
  logfile = './log/rnaseq-lib.' + today + '.log'
  log = logging.getLogger(__name__)
  log.setLevel(logging.INFO)
  log.addHandler(logging.StreamHandler())
  log.addHandler(logging.FileHandler(logfile))
  log.info('rnaseq-lib')
  log.info('Git version of rnaseq-lib:\n\t' + version)
  log.info('Runtime:\n\t' + today)
  log.info('Executed Python script: ')
  log.info('\t' + os.getcwd() + '/' +  __main__.__file__)
  log.info('Python:')
  temp = sys.version
  temp = temp.replace('\n','\n\t')
  log.info('\t' + temp)
  log.info('Packages:')
  log.info('\tPandas (v' + str(pandas.__version__) + ')')
  try:
    log.info('\tNumpy (v' + str(np.__version__) + ')')
  except:
    log.info('\tNumpy (not loaded)')
  try:
    log.info('\tmatplotlib (v' + str(matplotlib.__version__) + ')')
  except:
    log.info('\tmatplotlib (not loaded)')
  try:
    log.info('\tscipy (v' + str(scipy.__version__) + ')')
  except:
    log.info('\tscipy (not loaded)')
  log.info('')
  return log


def unit_tests():
  sample_xml = load_xml_summary()

  unittests.test_sample_key(SAMPLE_INFO, sample_xml)

def list_varscan_result_files():
  """
  List the varscan result files under 10_11_2014-varscan
  """

  files = []
  mutation_result_files = glob.glob('../10_11_2014-varscan/HL*')
  for file in mutation_result_files:
    if os.path.isdir(file):
      snp_file = glob.glob(file + '/snps.txt')
      files.append(snp_file[0])
  return files

def differential_expression_limma_voom(log, folder, groups, comparisons, max_zero_count_fraction=0.5):
  """
  Create the R script to perform differential expression comparing given sets of genes

  uses limma voom

  folder: folder to save results to
  samples: lists of lists of samples to compare against each other
           Should contian at least two lists
  groupnames: list of the group names
  max_zero_count_fraction: remove genes with more than this amount of samples
                           with zero read counts
  """

  import subprocess

  if not os.path.exists(folder):
    os.mkdir(folder)

  log.info('Performing differential expression on ' + str(len(comparisons)) + ' comparisons in ' + folder + '. Comparisons:')

  rscript = './' + folder + '/' + folder + '.R'
  rlog = './' + folder + '/' + folder + '.log'
  filehandle = open(rscript, 'w')

  filehandle.write('#\n')
  filehandle.write('# Comparison ' + folder + '\n')
  filehandle.write('# This code was autogenerated\n')
  filehandle.write('#\n')
  filehandle.write('plot_and_save_results <- function(results, name) {\n')
  filehandle.write('\t#\n')
  filehandle.write('\t# This function just plots some diagnostic plots and saves the results\n')
  filehandle.write('\t#\n')
  filehandle.write('\twrite.csv(results,paste(\"./' + folder + '/\",name,\"-results.csv\",sep=\"\"))\n')

  filehandle.write('\tpdf(paste(\"./' + folder + '/\",name,\"-pvalue-hist.pdf\",sep=\"\"))\n')
  filehandle.write('\thist(results$P.Value, breaks=20, col=\"grey\")\n')
  filehandle.write('\tdev.off()\n')

  filehandle.write('\tpdf(paste(\"./' + folder + '/\",name,\"-padj-hist.pdf\",sep=\"\"))\n')
  filehandle.write('\thist(results$adj.P.Val, breaks=20, col=\"grey\")\n')
  filehandle.write('\tdev.off()\n')

  filehandle.write('}\n\n')

  filehandle.write('library(\"limma\")\n')
  filehandle.write('library(\"edgeR\")\n')
  filehandle.write('library(\"RColorBrewer\")\n')
  filehandle.write('library(\"gplots\")\n')
  filehandle.write('library(\"genefilter\")\n\n')
  filehandle.write('setwd(\"' + os.getcwd() + '\")\n')

  filehandle.write('data <- read.csv(\"../9_30_2014-copydata/count-matrix.csv\",row.names=1, check.names=FALSE)\n')
  filehandle.write('y <- DGEList(counts=data)\n')
  filehandle.write('isexpr <- rowSums(cpm(y)>1) >= 32\n')
  filehandle.write('y <- y[isexpr,]\n')
  filehandle.write('y$samples$lib.size <- colSums(y$counts)\n')
  filehandle.write('y <- calcNormFactors(y)\n')

  #print the sample names

  count = 0
  for g in groups.keys():
    filehandle.write('group' + g + ' <- c(')
    for i in groups[g]:
      if i != groups[g][-1]:
        filehandle.write('\"' + i + '\",')
      else:
        filehandle.write('\"' + i + '\")\n')

    count += 1

  #write the DE analysis script for each comparison

  filehandle.write('\n\n')
  for c in comparisons:
    log.info('\t' + c[0] + '_vs_' + c[1])

    filehandle.write('y2 = y[,c(group' + c[0] + ', group' + c[1] + ')]\n')
    filehandle.write('samples <- as.factor(c(rep("1",length(group' + c[0] + ')),rep("2",length(group' + c[1] + '))))\n')
    filehandle.write('design <- model.matrix(~samples) \n')
    filehandle.write('v <- voom(y2,design) \n')
    filehandle.write('fit <- lmFit(v,design)\n')
    filehandle.write('fit <- eBayes(fit)\n')
    filehandle.write('options(digits=3)\n')
    filehandle.write('res = topTable(fit,coef=2,n=Inf,sort="p")\n')
    filehandle.write('plot_and_save_results(res, \"' + folder + '-' + c[0] + '_vs_' + c[1] + '\")\n\n')

  filehandle.close()

  cmd = ['Rscript',rscript]
  log.info('Running: ' + ' '.join(cmd))

  subprocess.call(cmd, stdout=open(rlog, "w"), stderr=open("/dev/null"))

def differential_expression(folder, groups, comparisons, max_zero_count_fraction=0.5):
  """
  Create the R script to perform differential expression comparing given sets of genes

  folder: folder to save results to
  samples: lists of lists of samples to compare against each other
           Should contian at least two lists
  groupnames: list of the group names
  max_zero_count_fraction: remove genes with more than this amount of samples
                           with zero read counts
  """

  import subprocess

  if not os.path.exists(folder):
    os.mkdir(folder)

  log.info('Performing differential expression on ' + str(len(comparisons)) + ' comparisons in ' + folder + '. Comparisons:')

  rscript = './' + folder + '/' + folder + '.R'
  rlog = './' + folder + '/' + folder + '.log'
  filehandle = open(rscript, 'w')

  filehandle.write('#\n')
  filehandle.write('# Comparison ' + folder + '\n')
  filehandle.write('# This code was autogenerated\n')
  filehandle.write('#\n')
  filehandle.write('plot_and_save_results <- function(results, name) {\n')
  filehandle.write('\t#\n')
  filehandle.write('\t# This function just plots some diagnostic plots and saves the results\n')
  filehandle.write('\t#\n')
  filehandle.write('\tdd = results[with(results,order(padj)),]\n')
  filehandle.write('\twrite.csv(dd,paste(\"./' + folder + '/\",name,\"-results.csv\",sep=\"\"))\n')

  filehandle.write('\tpdf(paste(\"./' + folder + '/\",name,\"-pvalue-hist.pdf\",sep=\"\"))\n')
  filehandle.write('\thist(res$pvalue, breaks=20, col=\"grey\")\n')
  filehandle.write('\tdev.off()\n')

  filehandle.write('\tpdf(paste(\"./' + folder + '/\",name,\"-padj-hist.pdf\",sep=\"\"))\n')
  filehandle.write('\thist(res$padj, breaks=20, col=\"grey\")\n')
  filehandle.write('\tdev.off()\n')

  filehandle.write('\tpdf(paste(\"./' + folder + '/\",name,\"-MA-plot.pdf\",sep=\"\"))\n')
  filehandle.write('\tplotMA(results,ylim=c(-1,1))\n')
  filehandle.write('\tdev.off()\n')

  filehandle.write('}\n\n')

  filehandle.write('library(\"DESeq2\")\n')
  filehandle.write('library(\"RColorBrewer\")\n')
  filehandle.write('library(\"gplots\")\n')
  filehandle.write('library(\"genefilter\")\n\n')

  filehandle.write('data <- read.csv(\"../9_30_2014-copydata/count-matrix.csv\",row.names=1, check.names=FALSE)\n')
  filehandle.write('temp <- apply(data,1,function(x) return(sum(x==0)/length(x)))\n')
  filehandle.write('data <- data[temp<=' + str(max_zero_count_fraction) + ',]\n\n')

  #print the sample names

  count = 0
  for g in groups.keys():
    filehandle.write('group' + g + ' <- c(')
    for i in groups[g]:
      if i != groups[g][-1]:
        filehandle.write('\"' + i + '\",')
      else:
        filehandle.write('\"' + i + '\")\n')

    count += 1

  #write the DE analysis script for each comparison

  filehandle.write('\n\n')
  for c in comparisons:
    log.info('\t' + c[0] + '_vs_' + c[1])

    filehandle.write('setwd(\"' + os.getcwd() + '\")\n')
    filehandle.write('data2 <- data[,c(group' + c[0] + ', group' + c[1] + ')]\n')
    filehandle.write('coldata <- data.frame(treatment=c(rep(\"1\",length(group' + c[0] + ')),rep(\"2\",length(group' + c[1] + '))), row.names=colnames(data2))\n')
    filehandle.write('countTable <- DESeqDataSetFromMatrix(countData=data2,colData=coldata,design=~treatment)\n')
    filehandle.write('result <- DESeq(countTable)\n')
    filehandle.write('res <- results(result)\n')
    filehandle.write('plot_and_save_results(res, \"' + folder + '-' + c[0] + '_vs_' + c[1] + '\")\n\n')

  filehandle.close()

  cmd = ['Rscript',rscript]
  log.info('Running: ' + ' '.join(cmd))

  subprocess.call(cmd, stdout=open(rlog, "w"), stderr=open("/dev/null"))

def get_fastq_files_and_sample():
  """
  Returns a dictionary where the keys are the sample names and
  the values are the paths to the fastq files

  """

  fastq_files = {}
  for i in rnaseq.DIRS:
    sample_xml = rnaseq.load_xml_summary()
    summary = sample_xml[i]
    sub_data = rnaseq.SAMPLE_INFO[rnaseq.SAMPLE_INFO['dir']==i]

    files = glob.glob(rnaseq.CLUSTER_SRCDIR + i + '/' + "*txt.gz")
    for j in files:
      sequence_number = os.path.split(j)[1].split('_')[1]
      temp = sub_data[sub_data['sequence #']==sequence_number]
      temp = temp['Sample #'].tolist()[0]
      if temp not in fastq_files:
        fastq_files[temp] = [j]
      else:
        fastq_files[temp].append(j)

  return fastq_files

def write_gct_file(filename, expression_data, descriptions=None):
  """
  Create the gct file. See the following for formatting
  http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats


  filename: filename to save to. must be .gct file type
  expression_data: pandas.DataFrame object (n by m). n being the number of
                   genes and m being the number of samples. Index must be
                   the gene names
  descriptions: vector of gene descriptions
  """

  assert filename.find('.gct')!=-1,"Pass a .gct filename"

  samples = expression_data.columns.tolist()
  genes = expression_data.index.tolist()
  if descriptions is None:
    descriptions = ['NA']*expression_data.shape[0]
  temp = ['NAME','DESCRIPTION'] + samples

  gct = pandas.DataFrame(columns=temp, index=expression_data.index)
  gct['NAME'] = genes
  gct['DESCRIPTION'] = descriptions
  for i in samples:
    gct[i] = expression_data[i]
  gct.to_csv('temp.txt',sep='\t',index=False)

  gct = open('temp.txt', 'r').read()
  fout = open(filename, 'w')
  fout.write('#1.2\n')
  fout.write(str(expression_data.shape[0]) + '\t' + str(len(samples)) + '\n')
  fout.write(gct)
  fout.close()

  os.system('rm temp.txt')

def perform_gsea_preRanked(folder, groups, comparisons, data, reference='msigdb.v4.0.symbols.gmt', limma=False):
  """
  Perform the GSEA PreRanked analysis
  """

  for c in comparisons:

    log.info('Performing GSEA-PreRanked on ' + c[0] + '_vs_' + c[1] + ' in ' + folder)

    rnk_file = './' + folder + '/' + c[0] + '_vs_' + c[1] + '.rnk'
    de_result_file_path = './' + folder + '/' + folder + '-' + c[0] + '_vs_' + c[1] + '-results.csv'

    write_rnk_file(
      filename=rnk_file,
      de_result_file_path=de_result_file_path,
      data=data,
      limma=limma
    )

    outdir = './' + folder + '/' + c[0] + '_vs_' + c[1]
    rpt_label = folder + '_' + c[0] + '_vs_' + c[1]
    _run_gseaPreRanked(outdir=outdir, rpt_label=rpt_label, rnk_file=rnk_file, reference=reference)

def write_rnk_file(filename, de_result_file_path, data, limma):
  """
  Creat the RNK file. See the following for formatting
  http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

  """

  assert filename.find('.rnk')!=-1,"Pass a .rnk filename"

  homologs, filtered_list = mouse_to_human(data.index.tolist())

  de_result_file = pandas.read_csv(de_result_file_path, header=0, index_col=0)

  log.info('\tTotal genes: ' + str(de_result_file.shape[0]))

  genes = data.index.tolist()
  de_genes = de_result_file.index.tolist()
  de_homologs = ['']*len(de_genes)
  for i in range(0, len(genes)):
    if genes[i] in de_genes:
      j = de_genes.index(genes[i])
      de_homologs[j] = homologs[i]

  de_result_file.index = de_homologs
  result_data = pandas.DataFrame(columns=['name','rank'])
  result_data['name'] = de_homologs
  if not limma:
      result_data['rank'] = de_result_file['log2FoldChange'].tolist()
  else:
      result_data['rank'] = de_result_file['logFC'].tolist()
  result_data = result_data[~pandas.isnull(result_data['name'])]
  result_data = result_data[~pandas.isnull(result_data['rank'])]

  log.info('\tFound ' + str(result_data.shape[0]) + ' homologs')

  result_data.to_csv(filename,index=False,header=False, sep='\t')

def _run_gseaPreRanked(outdir='', rpt_label='', rnk_file='', reference='msigdb.v4.0.symbols.gmt'):
  """
  Run the GSEA PreRanked analysis. First create the RNK file, then run the analysis
  """

  cmd = 'java -cp /Users/charlesmurphy/Downloads/gsea2-2.1.0.jar ' + \
        '-Xmx2048m xtools.gsea.GseaPreranked ' + \
        '-rnk ' + rnk_file + ' ' + \
        '-gmx ../../data/msigdb/' + reference + ' ' + \
        '-collapse false ' + \
        '-mode Max_probe ' + \
        '-norm meandiv ' + \
        '-nperm 1000 ' + \
        '-scoring_scheme classic ' + \
        '-rpt_label ' + reference + '-' + rpt_label + ' ' + \
        '-include_only_symbols true ' + \
        '-make_sets true ' + \
        '-plot_top_x 20 ' + \
        '-rnd_seed timestamp ' + \
        '-set_max 500 ' + \
        '-set_min 15 ' + \
        '-zip_report false ' + \
        '-out ' + outdir + ' ' + \
        '-gui false'
  log.info(cmd)
  os.system(cmd)

def run_gsea(prefix='', data=None, descriptions=None, class_labels=[]):
  """
  Create the .gct and .cls files then run gsea

  """

  if not os.path.exists(prefix):
    os.mkdir(prefix)

  write_gct_file(filename='./' + prefix + '/' + prefix + '.gct', expression_data=data, descriptions=descriptions)
  write_cls_file(filename='./' + prefix + '/' + prefix + '.cls', class_labels=class_labels)

  cmd = 'java -cp /Users/charlesmurphy/Downloads/gsea2-2.1.0.jar ' + \
        '-Xmx2048m xtools.gsea.Gsea ' + \
        '-res ' + prefix + '/' + prefix + '.gct' + \
        '-cls ' + prefix + '/' + prefix + '.cls#WT_versus_HET' + \
        '-gmx ../../data/msigdb/msigdb.v4.0.entrez.gmt' + \
        '-collapse false ' + \
        '-mode Max_probe ' + \
        '-norm meandiv ' + \
        '-nperm 20 ' + \
        '-permute phenotype ' + \
        '-rnd_type no_balance ' + \
        '-scoring_scheme weighted ' + \
        '-rpt_label my_analysis ' + \
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
        '-zip_report false ' + \
        '-out ./' + prefix + ' ' + \
        '-gui false'
  log.info(cmd)
  os.system(cmd)

"""
def write_GSEA_analysis(filename, output_dir):
  \"""
  Create the R script to run GSEA

  filename: R script filename
  \"""

  today = _get_time()

  fout = open(filename,'w')
  fout.write('# GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute\n')
  fout.write('#\n')
  fout.write('# R script to run GSEA Analysis \n')
  fout.write('#This script was created on: ' + today + '\n\n\n')

  fout.write('GSEA.program.location <- \"../10_27_14-gsea-r/GSEA.1.0.R\"\n')
  fout.write('source(GSEA.program.location, verbose=T, max.deparse.length=9999)\n')

  fout.write('GSEA(                                                                      # Input/Output Files :-------------------------------------------\n')
  fout.write(' input.ds =  \"/Users/charlesmurphy/Downloads/GSEA-P-R/Datasets/Gender.gct\",               # Input gene expression Affy dataset file in RES or GCT format\n')
  fout.write(' input.cls = \"/Users/charlesmurphy/Downloads/GSEA-P-R/Datasets/Gender.cls\",               # Input class vector (phenotype) file in CLS format\n')
  fout.write(' gs.db =     \"../../data/msigdb/msigdb.v4.0.entrez.gmt\",           # Gene set database in GMT format\n')
  fout.write(' output.directory      = \"' + output_dir + '\",            # Directory where to store output and results (default: "")\n')
  fout.write(' #  Program parameters :----------------------------------------------------------------------------------------------------------------------------\n')
  fout.write(' doc.string            = \"' + output_dir + '\",     # Documentation string used as a prefix to name result files (default: \"GSEA.analysis\")\n')
  fout.write(' non.interactive.run   = T,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)\n')
  fout.write(' reshuffling.type      = \"sample.labels\", # Type of permutation reshuffling: "sample.labels" or \"gene.labels\" (default: \"sample.labels\" \n')
  fout.write(' nperm                 = 1000,            # Number of random permutations (default: 1000)\n')
  fout.write(' weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)\n')
  fout.write(' nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)\n')
  fout.write(' fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)\n')
  fout.write(' fdr.q.val.threshold   = 0.25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)\n')
  fout.write(' topgs                 = 20,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)\n')
  fout.write(' adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)\n')
  fout.write(' gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)\n')
  fout.write(' gs.size.threshold.max = 500,             # Maximum size (in genes) for database gene sets to be considered (default: 500)\n')
  fout.write(' reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)\n')
  fout.write(' preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)\n')
  fout.write(' random.seed           = 111,             # Random number generator seed. (default: 123456)\n')
  fout.write(' perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)\n')
  fout.write(' fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)\n')
  fout.write(' replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)\n')
  fout.write(' save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)\n')
  fout.write(' OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)\n')
  fout.write(' use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)\n')
  fout.write(')\n')
  fout.close()
"""

def write_cls_file(filename, class_labels):
  """
  Create the cls file. See the following for formatting
  http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats

  filename: filename to save to. must be .cls file type
  class_labels: list of class labels for the expression data. Must be in the
                same order as the samples in the gct file
  """

  assert filename.find('.cls')!=-1,"Pass a .cls filename"
  classes = []
  for i in class_labels:
    if i not in classes:
      classes.append(i)

  fout = open(filename, 'w')
  fout.write(str(len(class_labels)) + '\t' + str(len(classes)) + '\t1\n')
  fout.write('#\t' + '\t'.join(classes) + '\n' + class_labels[0])
  for i in range(1, len(class_labels)):
    fout.write('\t' + class_labels[i])
  fout.write('\n')

  fout.close()


def mouse_to_human(mouse_list, idType='EntrezGene ID', recompute=False):
  """
  Convert a list of mouse gene symbols to human symbols using the
  MGI homolog database:
  ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt

  mouse_list: list of mouse identifiers
  idType: type of ID to use. Available: Symbol, EntrezGene ID, OMIM Gene ID, HGNC ID
  recompute: if True, call _comeute_mouse_human_homolog()
  """

  conversion_file = '../../lib/rnaseq-lib/data/mouse_to_human_homologs-' + idType + '.csv'
  if not os.path.exists(conversion_file) or recompute:
    _compute_mouse_human_homolog(mouse_list, idType='Symbol', conversion_file=conversion_file)
  try:
    save_data = pandas.read_csv(conversion_file, index_col=0, dtype={'human_homolog':str})
    homologs = save_data['human_homolog']
    filtered_mouse_list = save_data['mouse_ids']
  except:
    import pdb; pdb.set_trace()

  return homologs.tolist(), filtered_mouse_list.tolist()

def human_to_mouse(human_list, idType='EntrezGene ID', recompute=False):
  """
  Convert a list of mouse gene symbols to human symbols using the
  MGI homolog database:
  ftp://ftp.informatics.jax.org/pub/reports/HOM_MouseHumanSequence.rpt

  mouse_list: list of mouse identifiers
  idType: type of ID to use. Available: Symbol, EntrezGene ID, OMIM Gene ID, HGNC ID
  recompute: if True, call _comeute_mouse_human_homolog()
  """

  conversion_file = '../../lib/rnaseq-lib/data/human_to_mouse_homologs-' + idType + '.csv'
  #if not os.path.exists(conversion_file) or recompute:
  _compute_human_to_mouse_homolog(human_list, idType=idType, conversion_file=conversion_file)
  try:
    save_data = pandas.read_csv(conversion_file, index_col=0, dtype={'human_homolog':str})
    homologs = save_data['mouse_homolog']
    filtered_human_list = save_data['human_ids']
  except:
    import pdb; pdb.set_trace()

  return homologs.tolist(), filtered_human_list.tolist()


def load_xml_summary():
  sample_xml = {i:0 for i in DIRS}

  sample_xml['HuiLiu_2014_05_12'] = get_sample_key(LIBPATH + '/samples/HuiLiu_2014_05_12.xml')
  sample_xml['HuiLiu_2014_06_24'] = get_sample_key(LIBPATH + '/samples/HuiLiu_2014_06_24.xml')
  sample_xml['HuiLiu_2014_07_24'] = get_sample_key(LIBPATH + '/samples/HuiLiu_2014_07_24.xml')
  sample_xml['HuiLiu_2014_08_06'] = get_sample_key(LIBPATH + '/samples/HuiLiu_2014_08_06.xml')
  sample_xml['HuiLiu_2014_09_05'] = get_sample_key(LIBPATH + '/samples/HuiLiu_2014_09_05.xml')

  return sample_xml

def exit_handler():
  """
  Code that is executed upon program exit. It closes the loggin handlers
  """
  log.info("Done.")
  handlers = log.handlers[:]
  for handle in handlers:
      handle.close()
      log.removeHandler(handle)

def copy_files():
  """
  Copy the files from ole2001 directory to mine, keeping the same directory
  structure

  """

  log.info('In copy_files()')
  count = 0
  warnings = []

  for i in DIRS:

    summary_file = glob.glob(CLUSTER_SRCDIR + i + '/Summary.xml')[0]
    path_to_save = CLUSTER_WORKDIR + i + '/Summary.xml'
    shutil.copyfile(summary_file, path_to_save)

    subdirs = glob.glob(CLUSTER_SRCDIR + i + '/' + "*tophat")
    for j in subdirs:
      files = glob.glob(j + '/' + '*.bam.count')
      if len(files)==0:
        warnings.append(j)
        continue

      path_to_save = CLUSTER_WORKDIR + i
      if not os.path.exists(path_to_save):
        os.mkdir(path_to_save)

      path_to_save = CLUSTER_WORKDIR + i + '/' + os.path.split(j)[1]
      if not os.path.exists(path_to_save):
        os.mkdir(path_to_save)

      for file_to_copy in files:
        temp = path_to_save + '/' + os.path.split(file_to_copy)[1]
        if not os.path.exists(temp):
          shutil.copyfile(file_to_copy, temp)
          log.info('Copied to: ' + temp)
          count +=1
        else:
          log.info('Already exists: ' + temp)
  if len(warnings) != 0:
    log.info('WARNING: did not find .bam.count files for ' + str(len(warnings)) + ' samples')
    for w in warnings:
      log.info('\t' + str(w))
  log.info('Done, copied ' + str(count) + ' files.')

def get_sample_key(xml_file):
  """
  Load the Summary.xml file to get the lane-number, sample number key-value pairs
  """

  tree = ET.parse(xml_file)
  root = tree.getroot()
  temp_key = {}
  for lane in root.findall('Lane'):
    sampleId = lane.find('sampleId').text
    laneNumber = lane.find('laneNumber').text
    temp_key[sampleId] = laneNumber

  return temp_key

def list_bam_files():
  """
  Lists the bam files from ole2001's directory
  """

  log.info('In list_bam_files()')
  bam_files = []
  for i in DIRS:
    subdirs = glob.glob(CLUSTER_SRCDIR + i + '/' + "*tophat")
    for j in subdirs:
      files = glob.glob(j + '/' + 'accepted_hits.bam')
      bam_files += files
  return bam_files

def merge_into_matrix():
  """
  Merge the *.bam.count files into one count matrix
  """

  log.info('In merge_into_matrix()')
  sample_xml = load_xml_summary()

  #load the *.bam.count files

  data = {}
  for i in DIRS:
    subdirs = glob.glob(CLUSTER_WORKDIR + i + '/' + "*tophat")

    summary = sample_xml[i]
    sub_data = SAMPLE_INFO[SAMPLE_INFO['dir']==i]

    for j in subdirs:
      files = glob.glob(j + '/' + '*.bam.count')
      for file_to_load in files:

        sequence_number = os.path.split(j)[1].split('_')[1]
        temp = sub_data[sub_data['sequence #']==sequence_number]
        temp = temp['Sample ID'].tolist()
        assert len(temp)==1,"problem with in merge_into_matrix()"
        fin = open(file_to_load, 'r')
        data[temp[0]] = []
        for line in fin.readlines():
          line = line.rstrip()
          if line.find('__no_feature') != -1:
            break
          line = line.split('\t')
          data[temp[0]].append(line)
        fin.close()

  unittests.test_sample_load(data,log)

  #save to file

  names = data.keys()

  fout = open('count-matrix.csv', 'w')
  fout.write('gene')
  for i in names:
    fout.write(',' + i)
  fout.write('\n')
  for i in range(0, len(data[i])):
    fout.write(data[names[0]][i][0])
    for j in names:
      fout.write(',' + str(data[j][i][1]))
    fout.write('\n')
  fout.close()
  log.info('\tSave to count-matrix.csv (' + str(len(data[names[0]])) + ' rows, ' + str(len(names)) + ' columns)')
  log.info('\tSamples: ' + str(names))

def sort_by_expression_variance(data):
  """
  Sort the count matrix by variance
  """

  data2 = np.array(data)
  genes = data.index
  variance = data2.var(axis=1)
  sort_index = sorted(range(len(variance)), key=lambda k: variance[k], reverse=True)
  data2 = data2[sort_index]
  genes = genes[sort_index]

  return data2, genes

def remove_genes_zero_sample_counts(data, threshold=0.5, count_threshold=0):
  """
  Remove genes with zero counts in a given percent of samples
  """

  remaining_genes = []
  N = float(data.shape[1])

  for i in data.index:
    temp = (data.ix[i]<=0).sum() / N
    if temp <= threshold:
      remaining_genes.append(i)

  data = data.ix[remaining_genes]

  return data

def get_list_of_samples_by_number():
  """
  Returns list of samples by #. For example,

  [['1259-Control 4**','1259-BMN 2'],['1397-BMN-1','1397-BYL-2']]

  """

  samples = SAMPLE_INFO[SAMPLE_INFO['dir']!='na']

  parental = samples[samples['Parental']==1]
  non_parental = samples[samples['Parental']!=1]
  tumors = list(set(non_parental['Mouse']))

  samples = []
  for i in tumors:
    temp = non_parental[non_parental['Mouse']==i]
    temp = temp['Sample ID'].tolist()
    samples.append(temp)

  samples.append(parental['Sample ID'].tolist())

  return samples

def get_mouse_numbers(samples):
  """
  Given a list of samples names (e.g. [1513-Ctrl-2, 1397-BYL+BMN-1, ...])
  return a list of integers. Where each integer corresponds to mouse number
  or parental strain.

  """

  numbers = []
  index = {}
  count = 1
  for i in samples:
    if i.find('parental') == -1:
      temp = int(i.split('-')[0])
      if temp not in index:
        index[temp] = count
        count+=1
      numbers.append(index[temp])
    else:
      if 'parental' not in index:
        index['parental'] = count
        count+=1
      numbers.append(index['parental'])
  return numbers, index

def get_symbols_kegg_pathway(name):
    """
    Return the gene symbols for a particular KEGG pathway
    """

    fin = open('../../data/msigdb/c2.cp.kegg.v4.0.symbols.gmt','r')
    kegg = {}
    for line in fin.readlines():
        line = line.rstrip().split('\t')
        kegg[line[0]] = line[2:]
    fin.close()

    return kegg[name]

def get_expression_matrix_kegg_pathway(data,pathway=''):
    """
    Given a specific KEGG pathway, return the gene expression matrix containing
    the mouse gene homologs of the pathway

    data: the gene expression matrix
    pathway: the KEGG pathway
    """

    path = get_symbols_kegg_pathway(name=pathway)

    homologs, filtered_list = human_to_mouse(path, idType='Symbol', recompute=False)
    non_nan_homologs = []
    for i in homologs:
        if type(i)==str:
            non_nan_homologs.append(i)

    return data.ix[non_nan_homologs]

def jitter_plot_gene_expression(data,gene_name,prefix,samples,group_names,title='',rotation=0):
    """
    Plot jitter plot of gene across samples

    data: pandas.DataFrame object containing
          the n by m count matrix (n genes, m samples)
    gene_name: the name of the gene
    *samples: one or more list of samples of which to plot expression of the gene
    """

    assert gene_name in data.index, gene_name + " is not in the count matrix"
    assert len(samples)!=0, "Provide at least one list of sample names"


    all_samples = []
    for samp in samples:
      all_samples += samp
    number_sample_groups = len(samples)
    print all_samples

    x = []
    y = []
    for i in range(0, number_sample_groups):
        temp = data.ix[gene_name][samples[i]].tolist()
        y += temp
        for j in range(0, len(temp)):
            x.append(i+1 + (random.random()/20 - 0.025))

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_axes([0.1,0.2,0.8,0.7])
    barlist = ax.plot(x, y, 'o')

    # add some text for labels, title and axes ticks
    ax.set_ylabel('Expression',fontsize=22)
    ax.set_title('Expression of ' + gene_name,fontsize=22)
    ax.set_xticks(np.arange(number_sample_groups)+1)
    ax.set_xticklabels( group_names,fontsize=12,rotation=rotation)
    plt.title(title,fontsize=22)
    plt.savefig(prefix + '-' + gene_name + '.png')
    plt.clf()
    plt.close()

    log.info("jitter_plot_gene_expression()")
    log.info('\tGenes: ' + str(gene_name))
    log.info('\tSamples: ' + str(samples))

def plot_barplot_gene_expression(data, gene_name, prefix, samples):
  """
  Plot bar plot of gene across samples

  data: pandas.DataFrame object containing
        the n by m count matrix (n genes, m samples)
  gene_name: the name of the gene
  *samples: one or more list of samples of which to plot expression of the gene
  """

  assert gene_name in data.index, gene_name + " is not in the count matrix"
  assert len(samples)!=0, "Provide at least one list of sample names"


  all_samples = []
  for samp in samples:
    all_samples += samp
  number_sample_groups = len(samples)
  print all_samples

  exp = data.ix[gene_name][all_samples].tolist()
  ind = np.arange(len(all_samples))

  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.1,0.2,0.8,0.7])
  barlist = ax.bar(ind, exp, 0.25, color='r', edgecolor=None)

  # add some text for labels, title and axes ticks
  ax.set_ylabel('Expression')
  ax.set_title('Expression of ' + gene_name)
  ax.set_xticks(ind+0.25)
  ax.set_xticklabels( all_samples, rotation=90)

  #plot different colors for each sample group

  if number_sample_groups > 1:
    from pylab import get_cmap
    cm = get_cmap('Paired')

    colors = []
    for i in range(0, number_sample_groups):
      for j in range(0, len(samples[i])):
        colors.append(cm( i / float(number_sample_groups) ))

    for i in range(0, len(barlist)):
      barlist[i].set_color(colors[i])

    count = 0
    for ticklabel in ax.get_xticklabels():
      ticklabel.set_color(colors[count])
      count += 1

  plt.savefig(prefix + '-' + gene_name + '.png')
  plt.clf()
  plt.close()

  log.info("plot_barplot_gene_expression()")
  log.info('\tGenes: ' + str(gene_name))
  log.info('\tSamples: ' + str(samples))

def pca(data, color=False):
  """
  Perform PCA
  """

  samples = data.columns

  if color:
    mouse_number, mouse_number_key = get_mouse_numbers(samples)
    total = float(len(list(set(mouse_number))))
    cm = get_cmap('Paired')
    colors = []
    for i in range(0, len(mouse_number)):
      colors.append(cm(mouse_number[i] / total))
    #colors = np.array(colors)
  else:
    colors = ['b']*len(samples)

  pca = PCA(data)

  fig = plt.figure()
  fig.set_size_inches(12, 12)
  plt.scatter(pca.Y[:, 0], pca.Y[:, 1], c=colors, marker='o', edgecolor='none', s=40)
  plt.xlabel('PC1')
  plt.ylabel('PC2')
  plt.savefig('pca1-2.png')
  plt.clf()
  plt.close()

def pca_genes_as_features(data, color=False):
  """
  Perform PCA
  """
  samples = data.columns
  data, genes = sort_by_expression_variance(data)
  n = data.shape[1]


  if color:
    mouse_number, mouse_number_key = get_mouse_numbers(samples)
    total = float(len(list(set(mouse_number))))
    cm = get_cmap('Paired')
    colors = []
    for i in range(0, len(mouse_number)):
      colors.append(cm(mouse_number[i] / total))
    #colors = np.array(colors)
  else:
    colors = ['b']*len(samples)

  pca = PCA(data[0:n,].T)
  import pdb; pdb.set_trace()

  fig = plt.figure()
  fig.set_size_inches(12, 12)
  plt.scatter(pca.Y[:, 0], pca.Y[:, 1], c=colors, marker='o', edgecolor='none', s=40)
  plt.xlabel('PC1')
  plt.ylabel('PC2')
  plt.savefig('pca1-2-genesAsFeatures.png')
  plt.clf()
  plt.close()

def plot_correlation_heatmap(data, filename='heatmap-correlations.png', samples=[]):
  """
  plot a heatmap of the correlations between samples

  data: pandas.DataFrame object containing
        the n by m count matrix (n genes, m samples)
  filename: name to save image
  samples: a list containing specific samples to plot, default is to
           plot all samples

  """

  assert len(samples)!=0,"provide at least one list of samples"

  all_samples = []
  for samp in samples:
    all_samples += samp
  number_sample_groups = len(samples)

  data = data[all_samples]

  fig = plt.figure(figsize=(10,10))
  ax = fig.add_axes([0.2,0.2,0.7,0.7])
  heatmap = ax.pcolor(np.corrcoef(data.T), cmap=plt.cm.RdBu, alpha=0.8)
  plt.title('Correlation heatmap', size=16)
  ax.set_frame_on(False)
  ax.set_yticks(np.arange(data.T.shape[0]) + 0.5, minor=False)
  ax.set_yticklabels(data.T.index, minor=False, fontsize=8)

  ax.set_xticks(np.arange(data.T.shape[0]) + 0.5, minor=False)
  ax.set_xticklabels(data.T.index, minor=False, fontsize=8, rotation=90)

  #plot different colors for each sample group

  if number_sample_groups > 1:
    from pylab import get_cmap
    cm = get_cmap('Paired')

    colors = []
    for i in range(0, number_sample_groups):
      for j in range(0, len(samples[i])):
        colors.append(cm( i / float(number_sample_groups) ))

    count = 0
    for ticklabel in ax.get_xticklabels():
      ticklabel.set_color(colors[count])
      count += 1

    count = 0
    for ticklabel in ax.get_yticklabels():
      ticklabel.set_color(colors[count])
      count += 1


  ax.invert_yaxis()
  ax.grid(False)
  plt.colorbar(heatmap)
  plt.savefig(filename)
  plt.clf()
  plt.close()

  log.info("plot_correlation_heatmap()")
  log.info('\tSamples: ' + str(samples))

def load_data(file_path):
  """
  Load the matrix of RNA-Seq count data
  """
  log.info('Loading data...')

  data = pandas.read_csv(file_path, header=0, index_col=0)
  return data

def plot_hierarchical_clustering(data=[], filename = '', distance='euclidean', samples=[], fontsize=7):
  """

  data: pandas.DataFrame object containing
        the n by m count matrix (n genes, m samples)
  distance: distance measure to use
  """

  assert len(samples)!=0,"provide at least one list of samples"

  if filename == '':
    filename='hierarchical-clustering-' + distance + '.png'

  all_samples = []
  for samp in samples:
    all_samples += samp
  number_sample_groups = len(samples)

  if number_sample_groups > 1:
    from pylab import get_cmap
    cm = get_cmap('Paired')

    colors = []
    color_sample_map = {}
    for i in range(0, number_sample_groups):
      for j in range(0, len(samples[i])):
        colors.append(cm( i / float(number_sample_groups) ))


  data = data[all_samples]
  genes = data.index
  data = np.array(data).T

  fig = plt.figure(figsize=(20,16))

  data_dist = pdist(data, distance)
  ax1 = fig.add_axes([0.02,0.2,0.18,0.60])
  Y = linkage(data_dist, method='single')
  Z1 = dendrogram(Y, orientation='right',no_labels=True, color_threshold=0,distance_sort='descending')
  ax1.set_xticks([])

  data_dist = pdist(data, distance)
  ax2 = fig.add_axes([0.21,0.81,0.55,0.18])
  Y = linkage(data_dist, method='single')
  Z2 = dendrogram(Y, orientation='top',no_labels=True,color_threshold=0)
  ax2.set_xticks([])

  axmatrix = fig.add_axes([0.21,0.2,0.55,0.60])
  D = squareform(data_dist)
  idx1 = Z1['leaves']
  idx2 = Z2['leaves']
  D = D[idx1,:]
  D = D[:,idx2]
  labels1 = np.array(all_samples)[idx1]
  labels2 = np.array(all_samples)[idx2]
  im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.RdBu)
  axmatrix.xaxis.tick_bottom()
  axmatrix.set_xticks(np.arange(len(labels2)), minor=False)
  axmatrix.set_xticklabels(labels2, minor=False, fontsize=fontsize, rotation=90)
  axmatrix.yaxis.tick_right()
  axmatrix.set_yticks(np.arange(len(labels1)), minor=False)
  axmatrix.set_yticklabels(labels1, minor=False, fontsize=fontsize)
  #axmatrix.invert_yaxis()

  if number_sample_groups > 1:
    colors1 = np.array(colors)[idx1]
    colors2 = np.array(colors)[idx2]
    count = 0
    for ticklabel in axmatrix.get_xticklabels():
      ticklabel.set_color(colors2[count])
      count += 1

    count = 0
    for ticklabel in axmatrix.get_yticklabels():
      ticklabel.set_color(colors1[count])
      count += 1


  axcolor = fig.add_axes([0.81,0.81,0.02,0.18])
  plt.colorbar(im, cax=axcolor)
  plt.savefig(filename)
  plt.clf()
  plt.close()

  log.info("plot_hierarchical_clustering()")
  log.info('\tDistance measure: ' + distance)
  log.info('\tSamples: ' + str(samples))


def plot_expression_heatmap(data, filename='', samples=[], title='',xlabel_fontsize=10,ylabel_fontsize=10,ylabel_rotation=90):
  """
  plot heatmap of gene expression
  """

  all_samples = []
  for samp in samples:
    all_samples += samp
  number_sample_groups = len(samples)
  data = data[all_samples]

  if number_sample_groups > 1:
    from pylab import get_cmap
    cm = get_cmap('Paired')

    colors = []
    color_sample_map = {}
    for i in range(0, number_sample_groups):
      for j in range(0, len(samples[i])):
        colors.append(cm( i / float(number_sample_groups) ))

  fig = plt.figure(figsize=(14,10))
  ax = fig.add_axes([0.1,0.2,0.7,0.7])
  plt.text(0.5,1.1,title,fontsize=22,horizontalalignment='center',transform=ax.transAxes)
  heatmap = ax.pcolor(data.values, cmap=plt.cm.RdBu, alpha=0.8)
  fig = plt.gcf()
  ax.set_frame_on(False)
  ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor=False)
  ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor=False)
  ax.invert_yaxis()
  ax.grid(False)
  ax.set_yticklabels(data.index, minor=False, fontsize=ylabel_fontsize)
  ax.set_xticklabels(all_samples, minor=False, fontsize=xlabel_fontsize, rotation=ylabel_rotation)
  axcolor = fig.add_axes([0.81,0.71,0.02,0.18])
  plt.colorbar(heatmap, cax=axcolor)


  if number_sample_groups > 1:
    count = 0
    for ticklabel in ax.get_xticklabels():
      ticklabel.set_color(colors[count])
      count += 1

  plt.savefig(filename)
  plt.clf()

  log.info("plot_expression_heatmap()")
  log.info('\tNumber genes: ' + str(data.shape[0]))


def plot_top_exp_var(data, n = -1, filename=''):
  """
  plot top n genes with highest variance
  """

  if n == -1:
    n = data.shape[0]

  if filename=='':
    filename='top-' + str(n) + '-genes-by-variation.png'

  data2, genes = sort_by_expression_variance(data)

  fig, ax = plt.subplots()
  heatmap = ax.pcolor(data2[0:n,], cmap=plt.cm.RdBu, alpha=0.8)
  fig = plt.gcf()
  fig.set_size_inches(14, 12)
  ax.set_frame_on(False)
  ax.set_yticks(np.arange(data2[0:n,].shape[0]) + 0.5, minor=False)
  ax.set_xticks(np.arange(data2.shape[1]) + 0.5, minor=False)
  ax.invert_yaxis()
  ax.grid(False)
  ax.set_yticklabels(genes[0:n], minor=False, fontsize=8)
  ax.set_xticklabels(data.columns, minor=False, fontsize=8, rotation=90)
  plt.colorbar(heatmap)
  plt.savefig(filename)
  plt.clf()

  log.info("plot_top_exp_var()")
  log.info('\tNumber genes: ' + str(n))

def plot_top_exp_var_cluster(data, n = -1, filename='', distance='euclidean', yaxisLabels=True, samples=[], y_fontsize=8):
  """
  plot top n genes with highest variance
  """

  all_samples = []
  for samp in samples:
    all_samples += samp
  number_sample_groups = len(samples)
  data = data[all_samples]

  if n == -1:
    n = data.shape[0]

  if filename=='':
    filename='top-' + str(n) + '-genes-by-variation.png'

  data2, genes = sort_by_expression_variance(data)
  data2 = data2[0:n,]
  genes = genes[0:n]

  fig = plt.figure(figsize=(16,12))

  data_dist = pdist(data2.T, distance)
  ax1 = fig.add_axes([0.1,0.81,0.7,0.18])
  Y = linkage(data_dist, method='average')
  Z2 = dendrogram(Y, orientation='top',no_labels=True,color_threshold=0)
  ax1.set_xticks([])


  ax2 = fig.add_axes([0.1,0.2,0.7,0.60])
  idx2 = Z2['leaves']
  labels = np.array(all_samples)[idx2]

  data2 = data2[:,idx2]
  heatmap = ax2.matshow(data2, aspect='auto', origin='lower', cmap=plt.cm.RdBu)
  ax2.set_frame_on(False)

  if yaxisLabels:
    ax2.set_yticks(np.arange(data2.shape[0]), minor=False)
    ax2.set_yticklabels(genes, minor=False, fontsize=y_fontsize)
  else:
    ax2.set_yticks([])

  ax2.xaxis.tick_bottom()
  ax2.set_xticks(np.arange(len(labels)), minor=False)
  ax2.set_xticklabels(labels, minor=False, fontsize=8, rotation=90)
  ax2.invert_yaxis()

  ax2.grid(False)
  axcolor = fig.add_axes([0.81,0.81,0.02,0.18])
  plt.colorbar(heatmap, cax=axcolor)
  plt.savefig(filename)
  plt.clf()

  log.info("plot_top_exp_var()")
  log.info('\tNumber genes: ' + str(n))
  plt.close()

unit_tests()
log = setup_log(log)
atexit.register(exit_handler)

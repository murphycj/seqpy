
import os
import datetime
import glob
import shutil
import imp
import sys
import random

import pandas
import scipy
import numpy as np


class ConnectivityMap():
    def __init__(self):
        self.conversion = '../../data/probe.conversions.txt'

    def prep(self):

        conversion = {}
        cmap = pandas.read_table(self.conversion,sep='\t',header=None)
        for i in cmap.index.tolist():
            temp = cmap.ix[i][1]
            if temp not in conversion:
                conversion[temp] = [cmap.ix[i][0]]
            else:
                conversion[temp].append(cmap.ix[i][0])
        return conversion

    def human_to_cmap(data):
        cmap = prep()
        results = {}
        for gene in data.index.tolist():
            if gene in cmap:
                for probe in cmap[gene]:
                    results[probe] = data.ix[gene]['log2FoldChange'].mean()
        return results


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

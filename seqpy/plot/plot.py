"""
Plotting functions for high-throughput sequencing data

"""

import os
import git
import datetime
import glob
import sys
import random

import pandas
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.mlab import PCA
from pylab import get_cmap
import scipy
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import numpy as np

def jitter(data,gene_name,prefix,samples,group_names,title='',rotation=0):
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

def barplot(data, gene_name, prefix, samples):
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

def expression_heatmap(data=[], filename = '', distance='euclidean', samples=[], fontsize=7):
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

  plt.close()

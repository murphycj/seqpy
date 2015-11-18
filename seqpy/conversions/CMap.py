
import os
import glob
import sys
import random

import pandas
import scipy
import numpy as np

dirname, filename = os.path.split(__file__)
CONVERSION_FILE = os.path.join(dirname,'data','probe.conversions.txt')

def get_cmap_conversion_table(self):
    """
    Load the data to convert human gene symbols to probe ids

    """

    conversion = {}
    cmap = pandas.read_table(CONVERSION_FILE,sep='\t',header=None)
    for i in cmap.index.tolist():
        temp = cmap.ix[i][1]
        if temp not in conversion:
            conversion[temp] = [cmap.ix[i][0]]
        else:
            conversion[temp].append(cmap.ix[i][0])

    return conversion

def human_to_cmap(data):
    """
    Given a list of human gene symbols and some expression metric

    data: pandas dataframe where index are genes and has a single column
          representing the expression measure
    """

    cmap_conversion_table = get_cmap_conversion_table()

    expression_measure = data.columns
    results = {}
    for gene in data.index.tolist():
        if gene in cmap_conversion_table:
            for probe in cmap_conversion_table[gene]:
                results[probe] = data.ix[gene][expression_measure].mean()

    return results

"""
Parses CRISPResso output files
"""

import glob
import os
import sys
import pandas


class CrispRVariants(object):
    """docstring for ."""
    def __init__(self, outdir, sample=''):
        self.outdir = outdir
        self.sample = sample

        self.average = 0
        self.median = 0
        self.overall = 0
        self.readCount = 0

    def parse_efficiency(self):

        quntification_file = glob.glob(self.outdir + '/*mutationEfficiency.txt')

        if len(quntification_file) != 1:
            print 'Error in mutationEfficiency.txt file.'
            sys.exit()

        d = pandas.read_table(
            quntification_file[0],
            index_col=0,
            sep=' '
        )
        self.average = d.ix['Average']['x']
        self.median = d.ix['Median']['x']
        self.overall = d.ix['Overall']['x']
        self.readCount = d.ix['ReadCount']['x']

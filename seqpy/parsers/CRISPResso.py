"""
Parses CRISPResso output files
"""

import glob
import os
import sys
import re


class CRISPResso(object):
    """docstring for ."""
    def __init__(self, outdir, sample=''):
        self.outdir = outdir
        self.sample = sample

        self.unmodified = 0
        self.nhej = 0
        self.hdr = 0
        self.mixed = 0
        self.total = 0

    def parse_efficiency(self):

        quntification_file = os.path.join(self.outdir, 'Quantification_of_editing_frequency.txt')

        if not os.path.exists(quntification_file):
            print 'No Quantification_of_editing_frequency.txt file found!'
            sys.exit()

        fin = open(quntification_file, 'r')
        fin.next()
        line = fin.next().rstrip().lstrip()
        if re.findall('Unmodified', line):
            self.unmodified = int(re.findall('\d+', line)[0])
        else:
            print 'Something wrong ' + line
            sys.exti()

        line = fin.next().rstrip().lstrip()
        if re.findall('NHEJ', line):
            self.nhej = int(re.findall('\d+', line)[0])
        else:
            print 'Something wrong ' + line
            sys.exti()

        line = fin.next().rstrip().lstrip()
        if re.findall('HDR', line):
            self.hdr = int(re.findall('\d+', line)[0])
        else:
            print 'Something wrong ' + line
            sys.exti()

        line = fin.next().rstrip().lstrip()
        if re.findall('Mixed HDR-NHEJ', line):
            self.mixed = int(re.findall('\d+', line)[0])
        else:
            print 'Something wrong ' + line
            sys.exti()

        fin.next()

        line = fin.next().rstrip().lstrip()
        if re.findall('Total Aligned', line):
            self.total = int(re.findall('\d+', line)[0])
        else:
            print 'Something wrong ' + line
            sys.exti()

"""
Parses picardtools output files
"""

import re


class CollectInsertSizeMetrics(object):
    """docstring for ."""
    def __init__(self, filename, sample=''):
        self.filename = filename
        self.metrics = {}
        self._parse()

    def _parse(self):
        fin = open(self.filename, 'r')

        while not re.findall('^## METRICS CLASS', fin.next()):
            pass

        header = fin.next().rstrip().split('\t')
        data = fin.next().rstrip().split('\t')

        for i, j in zip(header, data):
            self.metrics[i] = j

        fin.next()
        fin.next()
        histogram_header = fin.next().rstrip().split('\t')
        histogram_data = []

        for line in fin:
            histogram_data.append(line.rstrip().split('\t'))

        fin.close()

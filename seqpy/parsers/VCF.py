"""
Tools to parse stuff from VCF files
"""

import vcf
import pkg_resources
import pickle
from sys import exit
from seqpy.parsers.SNPEff import SnpEffInfo


class SummarizeVCFAnnotation(object):
    """

    """

    def __init__(self, filename, sample='', annotation=None):
        self.filename = filename
        self.data = None
        self.samples = None
        if annotation is None:
            print 'ERROR! Provide annotation when initializing SummarizeVCFAnnotation!'
            exit()
        self.annotation = annotation

    def summarize(self,samples=None):
        vcf_in = vcf.Reader(open(self.filename, 'r'))
        if samples is not None:
            self.data = {i: [] for i in samples}
        else:
            self.data = {i: [] for i in vcf_in.samples}

        for variant in vcf_in:
            variant_info = SnpEffInfo(variant.INFO)

            if variant_info.has_ann():
                anns = []
                for ann in variant_info.ann:
                    anns += ann.annotation

                mutation_classification = self.annotation.classify(anns)

            else:
                if len(variant.REF)>1 or len(variant.ALT[0])>1:
                    mutation_classification = 'indel'
                else:
                    mutation_classification = 'synonymous'

                print 'ERROR! Variant at chromosome %s, position %s missing INFO!' % (variant.CHROM,variant.POS)
                print 'Marking its annotation as %s.' % mutation_classification

            for s in variant.samples:
                if s.called and s.sample in self.data:
                    self.data[s.sample] += [mutation_classification]

        self.samples = self.data.keys()

    def save_summary(self, outfile):
        fout = open(outfile, 'w')
        fout.write(
            '%s,%s,%s\n' %
            (
                'sample',
                'total',
                ','.join(self.annotation.classes)
            )
        )
        for i in self.samples:
            fout.write(
                '%s,%s,%s\n' %
                (
                    i,
                    str(len(self.data[i])),
                    ','.join(map(lambda x: str(self.data[i].count(x)),self.annotation.classes))
                )
            )

        fout.close()

class VAF:
    def __init__(self, sample=None, mutect=False, varscan=False, pindel=False):
        self.sample = sample
        self.has_AD = hasattr(self.sample.data, 'AD')
        self.has_FREQ = hasattr(self.sample.data, 'FREQ')
        self.has_AF = hasattr(self.sample.data, 'AF')
        self.mutect = mutect
        self.pindel = pindel
        self.varscan = varscan

        self.freq = 'NA'
        self.mutant = 'NA'
        self.reference = 'NA'

        self._parse()

    def _parse_mutect(self):
        tmp = self.sample.data.AD
        while None in tmp:
            tmp.remove(None)
        self.freq = self.sample.data.AF
        self.mutant = tmp[1]
        self.reference = sum(tmp) - self.mutant

    def _parse_varscan(self):
        tmp = self.sample.data.FREQ
        tmp = float(tmp.replace('%', ''))/100.
        self.freq = tmp
        self.mutant = self.sample.data.AD
        self.reference = self.sample.data.RD

    def _parse_pindel(self):
        self.freq = self.sample.data.AD[1]/float(sum(self.sample.data.AD))
        self.mutant = self.sample.data.AD[1]
        self.reference = self.sample.data.AD[0]

    def _parse(self):

        if self.pindel and self.mutect:
            if hasattr(self.sample.data, 'AF'):
                self._parse_mutect()
            else:
                self._parse_pindel()
        elif self.pindel and self.varscan:
            pass
        elif self.mutect:
            self._parse_mutect()
        elif self.varscan:
            self._parse_varscan()
        elif self.pindel:
            self._parse_pindel()

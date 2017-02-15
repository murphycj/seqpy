"""
Parsers for getting variant allele frequency from vcf file
"""


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
        self.reference = sum(tmp)

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

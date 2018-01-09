"""
Tools to parse stuff from VCF files
"""

class PileupVcfAllele(object):

    def __init__(self, vcf_variant):

        alternative_alleles = []

        if len(vcf_variant.REF)==1:

            for alternative in vcf_variant.ALT:
                if len(alternative)>1:
                    alternative = str(alternative)[1::]

                    alternative = '+' + str(len(alternative)) + alternative

                    alternative_alleles.append(alternative)
                else:
                    alternative_alleles.append(str(alternative))
        else:
            reference = str(vcf_variant.REF)[1::]

            reference = '-' + str(len(reference)) + reference

            alternative_alleles.append(reference)

        self.alternative_alleles = alternative_alleles


class Pileup(object):
    """
    store data on genomic variations (snps and small indels)
    """

    def __init__(self, pileup_data):
        pileup_data = pileup_data.split('\t')

        self.seq_id = pileup_data[0]
        self.position = int(pileup_data[1])
        self.reference_base = pileup_data[2]
        self.coverage = int(pileup_data[3])

        #read bases and qualities are not printed if coverage is 0

        if self.coverage==0:
            self.read_bases = ''
            self.read_qualities = ''
        else:
            if len(pileup_data)!=6:
                print(pileup_data)
                sys.exit('Did not find 6 columns in pileup data!')
            self.read_bases = pileup_data[4]
            self.read_qualities = pileup_data[5]

    def base_count(self, allele):
        count = self.read_bases.count(allele.lower()) + self.read_bases.count(allele.upper())
        if (allele.upper() == self.reference_base) or (allele.lower() == self.reference_base):
            count += self.read_bases.count('.')
            count += self.read_bases.count(',')
        return count

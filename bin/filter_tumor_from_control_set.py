"""
filter_tumor_from_control_set.py

Given a VCF file containing tumor and normal (control) samples and their
corresponding set of BAM files, remove mutations from the "somatic" call set
using the given parameters. Tailored for RNA-seq data

"""

import vcf
import multiprocessing
import sys
import pandas
import os
import argparse


class Pileup:
    """
    store data on genomic variations (snps and small indels)
    """

    def __init__(self, filename, pileup_data):

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
                print pileup_data
                print filename
                sys.exit('Did not find 6 columns in pileup data!')
            self.read_bases = pileup_data[4]
            self.read_qualities = pileup_data[5]

    def base_count(self, base):
        count = self.read_bases.count(base.lower()) + self.read_bases.count(base.upper())
        if (base.upper() == self.reference_base) or (base.lower() == self.reference_base):
            count += self.read_bases.count('.')
            count += self.read_bases.count(',')
        return count

def parse_vcf(vcf_file):

    mutant_sites = {}

    vcf_in = vcf.Reader(open(vcf_file,'r'))
    for v in vcf_in:

        mutant_sites[str(v.CHROM) + '-' + str(v.POS)] = map(lambda x: x.sequence,v.ALT)
    return mutant_sites

def get_filter_site_set(args,pileup_files,mutant_sites):

    normal_coverage_count = {i:0 for i in mutant_sites.keys()}
    normal_mutant_support_count = {i:0 for i in mutant_sites.keys()}

    for pileup_file, sample in pileup_files:
        for line in open(pileup_file,'r').readlines():
            pileup = Pileup(pileup_file,line.rstrip().split('\t'))

            site = str(pileup.seq_id) + '-' + str(pileup.position)

            #assert site in mutant_sites, import pdb; pdb.set_trace()
            if site not in mutant_sites:
                import pdb; pdb.set_trace()


            #count number occurances of mutant base
            #loop in case there are multiple alternative alleles
            coverage_count = 0
            mutant_support_count = 0
            for b in mutant_sites[site]:
                count = pileup.base_count(b)

                #does the site have enough coverage?

                if pileup.coverage >= args.minCoverage:
                    coverage_count = 1

                #how many bases support the mutation?

                if count >= args.maxSupporting:
                    mutant_support_count = 1

            normal_coverage_count[site] += coverage_count
            normal_mutant_support_count[site] += mutant_support_count

    return normal_coverage_count, normal_mutant_support_count

def main(args):

    mutant_sites = parse_vcf(args.vcf)

    pileup_files = zip(args.pileups,args.samples)

    normal_coverage_count, normal_mutant_support_count = get_filter_site_set(args,pileup_files,mutant_sites)

    #filter vcf

    vcf_in = vcf.Reader(open(args.vcf,'r'))
    vcf_out = vcf.Writer(open(args.out,'w'), vcf_in)

    for variant in vcf_in:
        site = str(variant.CHROM) + '-' + str(variant.POS)
        if normal_coverage_count[site]>=args.minNormalExpressed and normal_mutant_support_count[site]==0:
            vcf_out.write_record(variant)


parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,required=True)
parser.add_argument('--pileups',type=str,required=True,help='Pileup files for control samples',nargs='+')
parser.add_argument('--samples',type=str,required=True,help='Sample names (same order as pileups)',nargs='+')
parser.add_argument(
    '--maxSupporting',
    type=int,
    default=2,
    required=True,
    help='Maximum number of supporting mutant bases that can be observed in any normal samples.'
)
parser.add_argument(
    '--minCoverage',
    type=int,
    default=6,
    required=True,
    help='Minimum covered in normal sample sufficient for mutation calling.'
)
parser.add_argument(
    '--minNormalExpressed',
    type=int,
    required=True,
    default=1,
    help='Minimum number of normal samples with sufficient reads covering mutation site to consider whether the site is somatic in the tumors'
)
parser.add_argument('--out',type=str,help='Output file name.',required=True)
args = parser.parse_args()

main(args=args)

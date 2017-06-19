"""
filter_tumor_from_control_set.py

Given a VCF file containing tumor and normal (control) samples and their
corresponding set of BAM files, remove mutations from the "somatic" call set
using the given parameters. Tailored for RNA-seq data

"""

import vcf
import sys
import pandas
import os
import argparse
from seqpy.parsers import Mpileup

def parse_vcf(vcf_file):

    mutant_sites = {}

    vcf_in = vcf.Reader(open(vcf_file,'r'))
    for v in vcf_in:
        allele = Mpileup.PileupVcfAllele(vcf_variant=v)
        mutant_sites[str(v.CHROM) + '-' + str(v.POS)] = allele.alternative_alleles
    return mutant_sites

def get_filter_site_set(args,pileup_files,mutant_sites,coverage,support):

    for pileup_file, sample in pileup_files:
        for line in open(pileup_file,'r').readlines():
            pileup = Mpileup.Pileup(line.rstrip())

            site = str(pileup.seq_id) + '-' + str(pileup.position)


            coverage.ix[site][sample]=pileup.coverage
            mutant_support_count=0
            for b in mutant_sites[site]:
                mutant_support_count += pileup.base_count(b)

            support.ix[site][sample]=mutant_support_count


    support = support.fillna(0)
    coverage = coverage.fillna(0)

    coverage.to_csv('coverage.csv')
    support.to_csv('support.csv')

def main(args):

    mutant_sites = parse_vcf(args.vcf)

    coverage = pandas.DataFrame(index=mutant_sites.keys(),columns=args.samples)
    support = pandas.DataFrame(index=mutant_sites.keys(),columns=args.samples)

    pileup_files = zip(args.pileups,args.samples)

    get_filter_site_set(args,pileup_files,mutant_sites,coverage,support)


parser = argparse.ArgumentParser(description='')
parser.add_argument('--vcf',type=str,required=False,default=None)
parser.add_argument('--mutation',type=str,required=False,help="",default=None)
parser.add_argument('--pileups',type=str,required=True,help='Pileup files for control samples',nargs='+')
parser.add_argument('--samples',type=str,required=True,help='Sample names (same order as pileups)',nargs='+')
args = parser.parse_args()

main(args=args)

"""
Tools to parse stuff from VCF files
"""

import vcf
from sys import exit
import pkg_resources
import pickle

class _ANN:
    def __init__(self, ann):
        ann = ann.split('|')
        assert len(ann) == 16, 'not enough ANN fields'

        self.allele = ann[0]
        self.annotation = ann[1].split('&')
        self.putative_impact = ann[2]
        self.gene_name = ann[3]
        self.gene_id = ann[4]
        self.feature_type = ann[5]
        self.feature_id = ann[6]
        self.transcript_biotype = ann[7]
        self.rank = ann[8]
        self.basepair_change = ann[9]
        self.aminoacid_change = ann[10]
        self.cdna_position = ann[11]
        self.cds_position = ann[12]
        self.protein_position = ann[13]
        self.distance_to_feature = ann[14]
        self.errors = ann[15]

class SnpEffAnnotation:
    """
    Contains a representation of the sequence ontology that the SnpEff
    variant annotation is based off on. It also will classify a given
    variant's annotation, while properling handling conflicts
    """

    def __init__(self,
            classes=[
                    'missense','synonymous','indel','splicesite','frameshift',
                    'truncating','nonsense','other'],
            rankings=[
                ['nonsense','missense','truncating','frameshift','splicesite','indel','synonymous','UTR_variant','other'],
            ]):
        """
        classes:
            the classes of mutations to classify the mutations into

        rules:
            the ranking rules for classes in the case that a mutatoion is
            classified into multiple classes. The classifications on the
            left take higher priority.
        """

        self.classes = classes
        self.rankings = rankings
        self.ontology = pickle.load(
            open(
                pkg_resources.resource_filename(
                    'seqpy',
                    'data/070717.ontology.pickle'
                ),
                'rb'
            )
        )
        self.classification = {

            # missense

            'nonsynonymous_variant':['missense'],
            'initiator_codon_variant':['missense'],
            'terminator_codon_variant':['missense'],
            'start_lost':['missense'],

            # synonymous

            'intergenic_variant':['synonymous'],
            'intron_variant':['synonymous'],
            'non_coding_transcript_variant':['synonymous'],
            'synonymous_variant':['synonymous'],
            'UTR_variant':['synonymous'],
            'incomplete_transcript_variant':['synonymous'],
            'intragenic_variant':['synonymous'],
            'transcript_secondary_structure_variant':['synonymous'],
            'NMD_transcript_variant':['synonymous'],
            'translational_product_structure_variant':['synonymous'],
            'translational_product_function_variant':['synonymous'],
            'regulatory_region_variant':['synonymous'],
            'sequence_length_variant':['synonymous'],
            'start_retained_variant':['synonymous'],
            'intergenic_region':['synonymous'],

            # indel

            'inframe_variant':['indel'],
            'complex_transcript_variant':['indel'],

            # nonsense

            'stop_gained':['nonsense','missense','truncating'],

            # UTR

            'UTR_variant':['UTR_variant','synonymous'],

            # splice site

            'splicing_variant':['splicesite','synonymous'],

            # frameshift

            'frameshift_variant':['frameshift','indel'],

            # truncation

            'frameshift_truncation':['truncating','indel'],

            # other

            'structural_variant':['other','synonymous'],
            'functional_variant':['other','synonymous'],
            'sequence_feature':['other','synonymous'],
            'protein_protein_contact':['other','synonymous']
        }

    def classify(self,annotations):
        """
        For each variant annotation for a mutation, it will travel up the
        gene ontology tree until meeting one of the keys in
        self.classification and use the associated value as the
        classification
        """

        def get_next(aa,name):
            classes = []
            if name in self.classification and \
                    any([x in self.classes for x in self.classification[name]]):
                return self.classification[name]
            if name == 'sequence_variant':
                print 'WARN! Reached sequence_variant in classification for %s!' % aa
                return [name]
            try:
                for j in self.ontology[name]:
                    classes += get_next(aa,j)
            except KeyError:
                print 'ERROR! Could not find %s in ontology!' % name
                print 'Exiting...'
                exit()
            return classes

        classes = []

        for aa in annotations:
            classes += get_next(aa,aa)

        classes = list(set(classes))

        if len(classes)==1:
            return classes[0]
        elif len(classes)==0:
            print 'ERROR! No mutation class found!'
            print 'Exiting...'
            exit()

        # find the variant annotation classification that has the highest
        # rank based on self.ranking

        previous_count = len(classes)

        while len(classes)>1:
            pair_i = classes[0]
            pair_j = classes[1]

            for rank in self.rankings:
                if pair_i in rank and pair_j in rank:
                    if rank.index(pair_i) > rank.index(pair_j):
                        classes.remove(pair_i)
                    else:
                        classes.remove(pair_j)

            if previous_count == len(classes):
                print 'ERROR! Don\'t know how to rank %s and %s!' % (pair_i,pair_j)
                print 'Exiting...'
                exit()
            previous_count = len(classes)

        return classes[0]


class SnpEffInfo:

    def __init__(self,info):
        self.info=info
        self.ann = []

        if 'ANN' in self.info:
            for ann in self.info['ANN']:
                self.ann.append(_ANN(ann=ann))

    def has_effect(self,effect):
        for ann in self.ann:
            if effect in ann.annotation:
                return True
        return False

    def has_ann(self):
        if 'ANN' in self.info:
            return True
        else:
            return False

    def get_effects_by_gene(self,only_most_deleterious=True):
        r = {}
        for ann in self.ann:
            if ann.gene_name not in r:
                r[ann.gene_name] = ann.annotation
            else:
                r[ann.gene_name] += ann.annotation

        # report most

        for i,j in r.items():
            r[i] = list(set(j))

        if only_most_deleterious:
            for gene, effects in r.items():
                r[gene] = self.get_most_deleterious_effect()

        return r

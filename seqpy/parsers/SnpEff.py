EFFECTS = {
    'coding_sequence_variant':4,
    'chromosome':1,
    'inframe_insertion':2,
    'disruptive_inframe_insertion':2,
    'inframe_deletion':2,
    'disruptive_inframe_deletion':2,
    'downstream_gene_variant':4,
    'exon_variant':4,
    'exon_loss_variant':1,
    'frameshift_variant':1,
    'gene_variant':4,
    'intergenic_region':4,
    'conserved_intergenic_variant':4,
    'intragenic_variant':4,
    'intron_variant':4,
    'conserved_intron_variant':4,
    'miRNA':4,
    'missense_variant':2,
    'initiator_codon_variant':3,
    'stop_retained_variant':3,
    'rare_amino_acid_variant':1,
    'splice_acceptor_variant':1,
    'splice_donor_variant':1,
    'splice_region_variant':3,
    'stop_lost':1,
    'start_lost':1,
    'stop_gained':1,
    'synonymous_variant':3,
    'start_retained':3,
    'stop_retained_variant':3,
    'transcript_variant':4,
    'regulatory_region_variant':4,
    'upstream_gene_variant':4,
    '3_prime_UTR_variant':4,
    '3_prime_UTR_truncation + exon_loss':2,
    '5_prime_UTR_variant':4,
    '5_prime_UTR_truncation + exon_loss_variant':2,
    'sequence_feature + exon_loss_variant':2,
    'non_coding_exon_variant':4,
    '5_prime_UTR_premature_start_codon_gain_variant':3,
    'non_coding_transcript_exon_variant':4,
    'structural_interaction_variant':1,
    'sequence_feature':2,
    'exon_loss_variant':2,
    'conservative_inframe_insertion':3,
    'conservative_inframe_deletion':3,
    'non_coding_transcript_variant':4,
    'protein_protein_contact':1
}

class _ANN:
    def __init__(self,ann):
        ann = ann.split('|')
        assert len(ann)==16,'not enough ANN fields'

        self.allele = ann[0]
        self.annotation=ann[1].split('&')
        self.putative_impact=ann[2]
        self.gene_name=ann[3]
        self.gene_id=ann[4]
        self.feature_type=ann[5]
        self.feature_id=ann[6]
        self.transcript_biotype=ann[7]
        self.rank=ann[8]
        self.basepair_change=ann[9]
        self.aminoacid_change=ann[10]
        self.cdna_position=ann[11]
        self.cds_position=ann[12]
        self.protein_position=ann[13]
        self.distance_to_feature=ann[14]
        self.errors=ann[15]


class SnpEffInfo:

    def __init__(self,info):
        self.info=info
        self.ann=[]

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

    def get_most_deleterious_effect(self):

        max_effect = list()
        max_score = 5

        effects = map(lambda x: x.annotation,self.ann)
        effects = [item for sublist in effects for item in sublist]

        for effect in effects:
            score = EFFECTS[effect]
            if score < max_score:
                max_score=score
                max_effect=[effect]
            elif score == max_score:
                max_effect.append(effect)

        max_effect = list(set(max_effect))

        if len(max_effect)>1:
            print '!!!WARNING - there is more than one effect that are most deleterious %s' % str(max_effect)

        return max_effect

    def get_effects_by_gene(self,only_most_deleterious=True):
        r = {}
        for ann in self.ann:
            if ann.gene_name not in r:
                r[ann.gene_name] = ann.annotation
            else:
                r[ann.gene_name] += ann.annotation

        #report most

        for i,j in r.items():
            r[i]=list(set(j))

        if only_most_deleterious:
            for gene, effects in r.items():
                r[gene] = self.get_most_deleterious_effect()

        return r

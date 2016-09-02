
class _ANN:
    def __init__(self,ann):
        ann = ann.split('|')
        assert len(ann)==16,'not enough ANN fields'

        self.allele = ann[0]
        self.annotation=ann[1]
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
            if effect==ann.annotation:
                return True
        return False

    def has_ann(self):
        if 'ANN' in self.info:
            return True
        else:
            return False

    def get_most_deleterious_effect(self):
        for ann in self.ann:
            if ann.annotation=='missense_variant':
                return 'missense_variant'
            if ann.annotation=='splice_donor_variant&intron_variant':
                return 'splice_donor_variant&intron_variant'
            if ann.annotation=='frameshift_variant':
                return 'frameshift_variant'
            if ann.annotation=='stop_gained':
                return 'stop_gained'
            if ann.annotation=='start_lost':
                return 'start_lost'
            if ann.annotation=='disruptive_inframe_insertion' or ann.annotation=='inframe_deletion' or ann.annotation=='inframe_insertion':
                return 'inframe_indel'

        return 'synonymous'

    def get_effects_by_gene(self):
        r = {}
        for ann in self.ann:
            if ann.gene_name not in r:
                r[ann.gene_name] = [ann.annotation]
            else:
                r[ann.gene_name].append(ann.annotation)

        for i,j in r.items():
            r[i]=list(set(j))
            
        return r

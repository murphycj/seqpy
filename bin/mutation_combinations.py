import itertools
import pandas
import argparse

class Combinations():
    """
    Stores the genes and mutation combiations of interest.
    Loaded from a tab-separated file of the following formate

        GeneA GeneB GeneC
        MUT MUT MUT
        MUT WT MUT
        ...
        WT WT WT

    Or if all combinations are of interest

        GeneA GeneB GeneC


    """

    def __init__(self,filepath):
        self.filepath = filepath
        self.genes = []
        self.combinations = []
        self._load()

    def _load(self):
        """
        Parse the file
        """
        data = pandas.read_table(self.filepath,sep='\t')
        self.genes = data.columns.tolist()

        if data.shape[0]==0:
            combinations = ["".join(seq) for seq in itertools.product("TF",repeat=len(self.genes))]
            self.combinations = map(lambda x: map(lambda x: True if x=="T" else False,x), combinations)
        else:
            combinations = data.values.tolist()
            check = list(set(itertools.chain.from_iterable(combinations)))
            assert len(check)<=2, "only include MUT or WT in combiations file"
            self.combinations = map(lambda x: map(lambda x: True if x=="T" else False,x), combinations)


class MutationData():
    """
    Parent class for holding, parsing, and printing the mutation data
    """

    def __init__(self,combinations=Combinations,prefix='',samples_to_exclude=[]):
        self.combinations = combinations
        self.prefix = prefix
        self.samples_to_exclude = samples_to_exclude

    def print_combinations(self):
        fout = open(self.prefix + '.combinations.csv','w')
        fout.write(','.join(self.combinations.genes) + ',percent,n,total\n')
        total = data.shape[0]
        results = {
            'total':total,
            'tallies':{}
        }

        patient_stage_mutation = []

        for i in combinations:
            temp = map(lambda x: "MUT" if x==True else "WT",i)
            temp = ','.join(temp)

            #get patients in certain mutation combination

            data_sub = data_n[
                (data_n["APC"]!=i[0]) &
                (data_n["KRAS"]!=i[1]) &
                (data_n["BRAF"]!=i[2]) &
                (data_n["PIK3CA"]!=i[3]) &
                (data_n["TP53"]!=i[4])
                ]

            count = data_sub.shape[0]
            percent = float(data_sub.shape[0])/total
            results['tallies'][temp]={'count':count,'percent':percent}
            fout.write(temp + ',' + str(percent) + ',' + str(count) + ',' + str(total) + '\n')

        fout.close()

    def print_combinations_by_stage(self):
        #write patient stage information

        for i in data_sub.index.tolist():
            if i in clinical_stage:
                if clinical_stage[i]=='[Not Available]':
                    patient_stage_mutation.append([i,temp,"NA"])
                else:
                    patient_stage_mutation.append([i,temp,clinical_stage[i]])
            else:
                patient_stage_mutation.append([i,temp,"NA"])


    def count_by_stage(self):
        fout = open(filename,'w')
        stages = list(set(clinical_stage.values()))
        stages.sort()

        #count number of patietns in each stage
        total_stage_count = {i:0 for i in stages}
        for i in patient_mutation_data.index.tolist():
            if i in clinical_stage:
                total_stage_count[clinical_stage[i]]+=1
            else:
                total_stage_count["NA"]+=1

        fout.write(',,,,')
        for s in stages:
            fout.write(',' + s + ' (n=' + str(total_stage_count[s]) + '),')
        fout.write('\nAPC,KRAS,BRAF,PIK3CA,TP53')
        for s in stages:
            fout.write(',percent,count')
        fout.write('\n')
        for i in combinations:
            temp = map(lambda x: "MUT" if x==True else "WT",i)
            temp = ','.join(temp)
            fout.write(temp)
            for s in stages:
                count = 0
                for j in data:
                    if j[1]==temp and j[2]==s:
                        count+=1
                if count==0:
                    fout.write(',0,' + str(count))
                else:
                    fout.write(',' + str((float(count)/total_stage_count[s])) + ',' + str(count))
            fout.write('\n')
        fout.close()


class MutationDataCosmic(MutationData):
    """
    Class to parse mutation data from Cosmic
    """

    def __init__(self):
        MutationData.__init__(self)

class MutationDataICGC(MutationData):
    """
    Class to parse mutation data from ICGC
    """

    def __init__(self):
        MutationData.__init__(self)

    def parse_mutation_data(self,filename):
        data = pandas.read_table(filename,sep='\t',low_memory=False)

        samples = data['submitted_sample_id'].tolist()
        data['submitted_sample_id'] = map(lambda x: '-'.join(x.split('-')[0:3]),samples)

        indices_to_keep = []

        for i in data.index.tolist():
            if data.ix[i]['gene_affected'] in self.combinations.genes:
                indices_to_keep.append(i)

        data = data.ix[indices_to_keep]

        data = data[
            (data['consequence_type']=='stop_gained') |
            (data['consequence_type']=='disruptive_inframe_deletion') |
            (data['consequence_type']=='frameshift_variant') |
            (data['consequence_type']=='inframe_deletion') |
            (data['consequence_type']=='missense_variant') |
            (data['consequence_type']=='stop_lost') |
            (data['consequence_type']=='start_lost') |
            (data['consequence_type']=='disruptive_inframe_insertion')
        ]

        #filter out any specified samples

        filtered_indices = []
        for i in data.index.tolist():
            temp = data.ix[i]['submitted_sample_id']
            if temp not in samples_to_exclude:
                filtered_indices.append(i)
        data = data.ix[filtered_indices]

        #load patient stage information

        clinical = pandas.read_table(clinical_file,sep='\t')
        clinical_stage = dict(zip(clinical['icgc_donor_id'],clinical['donor_tumour_stage_at_diagnosis_supplemental']))

        #count mutation combinations

        samples = list(set(data['icgc_donor_id']))
        genes=['APC','KRAS','BRAF','PIK3CA','TP53']
        mutations = pandas.DataFrame(index=samples,columns=genes)
        mutations.values.fill(False)

        for i in mutations.index:
            temp = data[data['icgc_donor_id']==i]

            if 'ENSG00000134982' in temp['gene_affected'].tolist():
                mutations['APC'][i]=True

            if 'ENSG00000133703' in temp['gene_affected'].tolist():
                mutations['KRAS'][i]=True

            if 'ENSG00000157764' in temp['gene_affected'].tolist():
                mutations['BRAF'][i]=True

            if 'ENSG00000121879' in temp['gene_affected'].tolist():
                mutations['PIK3CA'][i]=True

            if 'ENSG00000141510' in temp['gene_affected'].tolist():
                mutations['TP53'][i]=True

        fout = open(prefix + '.combinations.csv','w')
        fout.write('APC,KRAS,BRAF,PIK3CA,TP53,percent,n,total\n')
        total = mutations.shape[0]
        results = {
            'total':total,
            'tallies':{}
        }

        patient_stage_mutation = []

        for i in combinations:
            temp = map(lambda x: "MUT" if x==True else "WT",i)
            temp = ','.join(temp)
            data_sub = mutations[
                (mutations["APC"]==i[0]) &
                (mutations["KRAS"]==i[1]) &
                (mutations["BRAF"]==i[2]) &
                (mutations["PIK3CA"]==i[3]) &
                (mutations["TP53"]==i[4])
                ]

            count = data_sub.shape[0]
            percent = float(count)/total
            results['tallies'][temp]={'count':count,'percent':percent}
            fout.write(temp + ',' + str(percent) + ',' + str(count) + ',' + str(total) + '\n')

            for i in data_sub.index.tolist():
                if i in clinical_stage:
                    if clinical_stage[i]=='[Not Available]':
                        patient_stage_mutation.append([i,temp,"NA"])
                    else:
                        patient_stage_mutation.append([i,temp,clinical_stage[i]])
                else:
                    patient_stage_mutation.append([i,temp,"NA"])

        fout.close()

        count_by_stage(patient_stage_mutation,prefix + '.combinations.by.stage.csv',clinical_stage,mutations)

        samples = list(set(data['submitted_sample_id']))

        return results,samples

class MutationDatacBIO(MutationData):
    """
    Class to parse mutation data from cBio
    """

    def __init__(self):
        MutationData.__init__(self)
        self.mutation_data = None
        self.clinical_data = None
        self.tumor_stage = None

    def parse_mutation_data(self,filename):
        data = pandas.read_table(filename,sep='\t')
        del data["COMMON"]
        temp = filter(lambda x: 'Unnamed' in x, data.columns.tolist())
        del data[temp[0]]
        data.index = data["GENE_ID"]
        del data["GENE_ID"]
        data = data.T
        data_n = data.isnull()

        #simplify the tcga namse (e.g. TCGA-A1-1021-C1 becomes TCGA-A1-1021)

        new_names = []
        for i in data_n.index.tolist():
            if i.find('TCGA')!=-1:
                new_names.append('-'.join(i.split('-')[0:3]))
            else:
                new_names.append(i)
        data_n.index = new_names

        self.mutation_data = data_n

    def parse_tumor_stage(self):

        self.clinical_data['ajcc_pathologic_tumor_stage'] = self.clinical_data['ajcc_pathologic_tumor_stage'].replace('[Not Available]','NA')
        self.clinical_data['ajcc_pathologic_tumor_stage'] = self.clinical_data['ajcc_pathologic_tumor_stage'].replace('[Discrepancy]','NA')
        self.clinical_data = self.clinical_data.drop([0,1])
        self.tumor_stage = dict(zip(clinical['bcr_patient_barcode'],self.clinical_data['ajcc_pathologic_tumor_stage']))

    def parse_clinical(self,filename):

        #load patient clinical information and get stage
        self.clinical_data = pandas.read_csv(filename,sep='\t')

        self.parse_tumor_stage()


class MutationDataTCGA(MutationData):
    """
    Class to parse mutation data from TCGA
    """

    def __init__(self):
        MutationData.__init__(self)




def compare_samples():
    ICGC_READ_US = pandas.read_table('simple_somatic_mutation.open.READ-US.tsv',sep='\t',low_memory=False)
    ICGC_READ_US_samples = list(set(ICGC_READ_US['submitted_sample_id']))
    ICGC_READ_US_samples = map(lambda x: '-'.join(x.split('-')[0:3]),ICGC_READ_US_samples)

    ICGC_COCA_CN = pandas.read_table('simple_somatic_mutation.open.COCA-CN.tsv',sep='\t',low_memory=False)
    ICGC_COCA_CN_samples = list(set(ICGC_COCA_CN['submitted_sample_id']))
    ICGC_COCA_CN_samples = map(lambda x: '-'.join(x.split('-')[0:3]),ICGC_COCA_CN_samples)

    ICGC_COAD_US = pandas.read_table('simple_somatic_mutation.open.COAD-US.tsv',sep='\t',low_memory=False)
    ICGC_COAD_US_samples = list(set(ICGC_COAD_US['submitted_sample_id']))
    ICGC_COAD_US_samples = map(lambda x: '-'.join(x.split('-')[0:3]),ICGC_COAD_US_samples)

    ICGC_samples = set(ICGC_READ_US_samples + ICGC_COCA_CN_samples + ICGC_COAD_US_samples)

    genentech, samples = main("coadread_genentech_coadread_genentech_mutations.txt","Genentech")
    genentech_samples = set(samples)
    mskcc, samples = main("coadread_mskcc_coadread_mskcc_mutations.txt","MSKCC")
    mskcc_samples = set(samples)
    tcga_provisional, samples = main('coadread_tcga_pub_coadread_tcga_pub_mutations.txt','tcga_provisional')
    tcga_samples = set(samples)

    import pdb; pdb.set_trace()


def cosmic_mutations_genome(filepath,prefix,samples_to_exclude):
    cosmic = pandas.read_csv(filepath)

    #exclude samples from pubmed id 24951259 because they are targeted sequencing, which did not include PI3KCA or BRAF

    samples_to_exclude += ['P' + str(x) for x in range(1,161)]

    #reformat any existing TCGA samples so to include only the first three items

    new_sample_names=[]
    for i in cosmic['Sample Name'].tolist():
        if i.find('TCGA')!=-1:
            temp = '-'.join(i.split('-')[0:3])
            new_sample_names.append(temp)
        else:
            new_sample_names.append(i)
    cosmic['Sample Name'] = new_sample_names

    #remove all samples previously found in TCGA and ICGC datasets

    for i in samples_to_exclude:
        cosmic = cosmic[cosmic['Sample Name']!=i]

    samples = list(set(cosmic['COSMIC Sample ID'].tolist()))
    genes=['APC','KRAS','BRAF','PIK3CA','TP53']
    mutations = pandas.DataFrame(index=samples,columns=genes)
    mutations.values.fill(False)

    #count mutation combinations

    for i in mutations.index:
        temp = cosmic[cosmic['COSMIC Sample ID']==i]
        temp = temp[~temp['COSMIC Mutation ID'].isnull()]
        if 'APC' in temp['Gene Name'].tolist():
            mutations['APC'][i]=True

        if 'KRAS' in temp['Gene Name'].tolist():
            mutations['KRAS'][i]=True

        if 'BRAF' in temp['Gene Name'].tolist():
            mutations['BRAF'][i]=True

        if 'PIK3CA' in temp['Gene Name'].tolist():
            mutations['PIK3CA'][i]=True

        if 'TP53' in temp['Gene Name'].tolist():
            mutations['TP53'][i]=True

    #save the data

    fout = open(prefix + '.combinations.csv','w')
    fout.write('APC,KRAS,BRAF,PIK3CA,TP53,percent,n,total\n')
    total = mutations.shape[0]
    results = {
        'total':total,
        'tallies':{}
    }
    for i in combinations:
        temp = map(lambda x: "MUT" if x==True else "WT",i)
        temp = ','.join(temp)
        data_sub = mutations[
            (mutations["APC"]==i[0]) &
            (mutations["KRAS"]==i[1]) &
            (mutations["BRAF"]==i[2]) &
            (mutations["PIK3CA"]==i[3]) &
            (mutations["TP53"]==i[4])
            ]
        count = data_sub.shape[0]
        percent = float(count)/total
        results['tallies'][temp]={'count':count,'percent':percent}
        fout.write(temp + ',' + str(percent) + ',' + str(count) + ',' + str(total) + '\n')
    fout.close()

    cosmic.to_csv(prefix + '.filtered.csv')
    return results

def cosmic_mutations_non_genome(filepath,prefix):
    cosmic = pandas.read_csv(filepath)
    cosmic = cosmic[cosmic['Whole gene screen']=='y']

    #reformat any existing TCGA samples so to include only the first three items
    new_sample_names=[]
    for i in cosmic['Sample Name'].tolist():
        if i.find('TCGA')!=-1:
            temp = '-'.join(i.split('-')[0:3])
            new_sample_names.append(temp)
        else:
            new_sample_names.append(i)
    cosmic['Sample Name'] = new_sample_names

    #get the list of all samples queried for all 5 genes

    samples_with_all_five = []
    indices = []
    g = cosmic.groupby("COSMIC Sample ID")
    for i in g:
        num_genes = len(set(i[1]['Gene Name'].tolist()))
        if num_genes==5:
            samples_with_all_five += i[1]['COSMIC Sample ID'].tolist()
            indices += i[1].index.tolist()
    samples_with_all_five = list(set(samples_with_all_five))
    indices.sort()

    genes=['APC','KRAS','BRAF','PIK3CA','TP53']
    mutations = pandas.DataFrame(index=samples_with_all_five,columns=genes)
    mutations.values.fill(False)

    #count mutation combinations

    num_genes = []
    g = cosmic.groupby("COSMIC Sample ID")
    for i in g:
        num_genes = len(set(i[1]['Gene Name'].tolist()))
        if num_genes==5:
            temp = i[1]
            temp = temp[~temp['COSMIC Mutation ID'].isnull()]
            sample = list(set(i[1]['COSMIC Sample ID'].tolist()))
            temp2 = []
            if 'APC' in temp['Gene Name'].tolist():
                mutations['APC'][sample]=True
                temp2.append(True)
            else:
                temp2.append(False)
            if 'KRAS' in temp['Gene Name'].tolist():
                mutations['KRAS'][sample]=True
                temp2.append(True)
            else:
                temp2.append(False)

            if 'BRAF' in temp['Gene Name'].tolist():
                mutations['BRAF'][sample]=True
                temp2.append(True)
            else:
                temp2.append(False)

            if 'PIK3CA' in temp['Gene Name'].tolist():
                mutations['PIK3CA'][sample]=True
                temp2.append(True)
            else:
                temp2.append(False)

            if 'TP53' in temp['Gene Name'].tolist():
                mutations['TP53'][sample]=True
                temp2.append(True)
            else:
                temp2.append(False)

            #if temp2==[False,False,False,False,False]:
            #    print sample

    #save the data

    fout = open(prefix + '.combinations.csv','w')
    fout.write('APC,KRAS,BRAF,PIK3CA,TP53,percent,n,total\n')
    total = mutations.shape[0]
    results = {
        'total':total,
        'tallies':{}
    }

    for i in combinations:
        temp = map(lambda x: "MUT" if x==True else "WT",i)
        temp = ','.join(temp)
        data_sub = mutations[
            (mutations["APC"]==i[0]) &
            (mutations["KRAS"]==i[1]) &
            (mutations["BRAF"]==i[2]) &
            (mutations["PIK3CA"]==i[3]) &
            (mutations["TP53"]==i[4])
            ]
        count = data_sub.shape[0]
        percent = float(count)/total
        results['tallies'][temp]={'count':count,'percent':percent}
        fout.write(temp + ',' + str(percent) + ',' + str(count) + ',' + str(total) + '\n')
    fout.close()

    cosmic.ix[indices].to_csv(prefix + '.filtered.csv')
    return results


def ICGC_main(filepath='',prefix='',samples_to_exclude=[],clinical_file=''):
    data = pandas.read_table(filepath,sep='\t',low_memory=False)

    samples = data['submitted_sample_id'].tolist()
    data['submitted_sample_id'] = map(lambda x: '-'.join(x.split('-')[0:3]),samples)

    data = data[
        (data['gene_affected']=='ENSG00000133703') | #KRAS
        (data['gene_affected']=='ENSG00000157764') | #BRAF
        (data['gene_affected']=='ENSG00000134982') | #APC
        (data['gene_affected']=='ENSG00000121879') | #PIK3CA
        (data['gene_affected']=='ENSG00000141510')
        ]

    data = data[
        (data['consequence_type']=='stop_gained') |
        (data['consequence_type']=='disruptive_inframe_deletion') |
        (data['consequence_type']=='frameshift_variant') |
        (data['consequence_type']=='inframe_deletion') |
        (data['consequence_type']=='missense_variant') |
        (data['consequence_type']=='stop_lost') |
        (data['consequence_type']=='start_lost') |
        (data['consequence_type']=='disruptive_inframe_insertion')
    ]

    #filter out any specified samples

    filtered_indices = []
    for i in data.index.tolist():
        temp = data.ix[i]['submitted_sample_id']
        if temp not in samples_to_exclude:
            filtered_indices.append(i)
    data = data.ix[filtered_indices]

    #load patient stage information

    clinical = pandas.read_table(clinical_file,sep='\t')
    clinical_stage = dict(zip(clinical['icgc_donor_id'],clinical['donor_tumour_stage_at_diagnosis_supplemental']))

    #count mutation combinations

    samples = list(set(data['icgc_donor_id']))
    genes=['APC','KRAS','BRAF','PIK3CA','TP53']
    mutations = pandas.DataFrame(index=samples,columns=genes)
    mutations.values.fill(False)

    for i in mutations.index:
        temp = data[data['icgc_donor_id']==i]

        if 'ENSG00000134982' in temp['gene_affected'].tolist():
            mutations['APC'][i]=True

        if 'ENSG00000133703' in temp['gene_affected'].tolist():
            mutations['KRAS'][i]=True

        if 'ENSG00000157764' in temp['gene_affected'].tolist():
            mutations['BRAF'][i]=True

        if 'ENSG00000121879' in temp['gene_affected'].tolist():
            mutations['PIK3CA'][i]=True

        if 'ENSG00000141510' in temp['gene_affected'].tolist():
            mutations['TP53'][i]=True

    fout = open(prefix + '.combinations.csv','w')
    fout.write('APC,KRAS,BRAF,PIK3CA,TP53,percent,n,total\n')
    total = mutations.shape[0]
    results = {
        'total':total,
        'tallies':{}
    }

    patient_stage_mutation = []

    for i in combinations:
        temp = map(lambda x: "MUT" if x==True else "WT",i)
        temp = ','.join(temp)
        data_sub = mutations[
            (mutations["APC"]==i[0]) &
            (mutations["KRAS"]==i[1]) &
            (mutations["BRAF"]==i[2]) &
            (mutations["PIK3CA"]==i[3]) &
            (mutations["TP53"]==i[4])
            ]

        count = data_sub.shape[0]
        percent = float(count)/total
        results['tallies'][temp]={'count':count,'percent':percent}
        fout.write(temp + ',' + str(percent) + ',' + str(count) + ',' + str(total) + '\n')

        for i in data_sub.index.tolist():
            if i in clinical_stage:
                if clinical_stage[i]=='[Not Available]':
                    patient_stage_mutation.append([i,temp,"NA"])
                else:
                    patient_stage_mutation.append([i,temp,clinical_stage[i]])
            else:
                patient_stage_mutation.append([i,temp,"NA"])

    fout.close()

    count_by_stage(patient_stage_mutation,prefix + '.combinations.by.stage.csv',clinical_stage,mutations)

    samples = list(set(data['submitted_sample_id']))

    return results,samples


#compare_samples()
samples_to_exclude = []
"""
genentech, samples = main("coadread_genentech_coadread_genentech_mutations.txt","4.15.16-Genentech")
samples_to_exclude += samples
print "Genentech"

mskcc, samples = main("coadread_mskcc_coadread_mskcc_mutations.txt","4.15.16-MSKCC")
samples_to_exclude += samples
print "MSKCC"

tcga, samples = main("coadread_tcga_coadread_tcga_mutations.txt","4.15.16-TCGA",COAD=True)
samples_to_exclude += samples
print "TCGA-COAD-USA"

"""

#tcga_provisional, samples = main('coadread_tcga_pub_coadread_tcga_pub_mutations.txt','tcga_provisional')
#samples_to_exclude += samples

ICGC_CN, samples = ICGC_main('simple_somatic_mutation.open.COCA-CN.tsv','4.15.16-COCA-CN',[],clinical_file='donor.COCA-CN.tsv')
samples_to_exclude += samples
print 'ICGC - CN'

ICGC_COAD, samples = ICGC_main('simple_somatic_mutation.open.COAD-US.tsv','4.15.16-COAD-US',samples_to_exclude,clinical_file='donor.COAD-US.tsv')
samples_to_exclude += samples
print 'ICGC - COAD'

ICGC_READ, samples = ICGC_main('simple_somatic_mutation.open.READ-US.tsv','4.15.16-READ-US',samples_to_exclude,clinical_file='donor.READ-US.tsv')
samples_to_exclude += samples
print 'ICGC - READ'

import pdb; pdb.set_trace()

COSMIC_genome = cosmic_mutations_genome('4.1.16-whole-genome-screen.csv','4.15.16-cosmic.whole-genome-screen',samples_to_exclude)

COSMIC_nonGenome = cosmic_mutations_non_genome('4.1.16-non-whole-genome-screen.csv','4.15.16-cosmic.non-whole-genome-screen')

#write combined results

fout = open('combined.combinations.csv','w')
fout.write('APC,KRAS,BRAF,PIK3CA,TP53,percent,n,total\n')
total = genentech['total'] + mskcc['total'] + tcga['total'] + ICGC_CN['total'] + ICGC_COAD['total'] + ICGC_READ['total'] + COSMIC_genome['total'] + COSMIC_nonGenome['total']
for i in combinations:
    temp = map(lambda x: "MUT" if x==True else "WT",i)
    temp = ','.join(temp)
    count = genentech['tallies'][temp]['count'] + mskcc['tallies'][temp]['count'] + tcga['tallies'][temp]['count'] + ICGC_CN['tallies'][temp]['count'] + ICGC_COAD['tallies'][temp]['count'] + ICGC_READ['tallies'][temp]['count'] + COSMIC_genome['tallies'][temp]['count'] + COSMIC_nonGenome['tallies'][temp]['count']
    percent = float(count)/total
    fout.write(temp + ',' + str(percent) + ',' + str(count) + ',' + str(total) + '\n')
fout.close()

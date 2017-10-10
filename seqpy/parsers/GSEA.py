import xlsxwriter
import pandas
import os
import numpy as np
import glob


class GSEA():
    """
    class to store and parse GSEA results

    """

    def __init__(self, directory='', name='', phenotype1='up', phenotype2='down'):
        """
        directory: the directory containing the full GSEA result files

        name: name to give the analysis, this will be the name on the
              tab in the excel file
        """

        self.directory = directory
        self.phenotype1 = phenotype1
        self.phenotype2 = phenotype2

    def _check_for_nas(self, table):
        """
        Sometimes GSEA shows that significant pathways have NA for p-value
        """

        if sum(table['NOM p-val'].isnull()) > 0:
            indices = table.index.tolist()

            print 'NAs found for NOM p-val for significant pathways in:'
            print '\t\t' + self.directory

            for i in range(0, len(table.index)):
                if pandas.isnull(table.ix[indices[i]]['NOM p-val']):
                    were_there_nas = True
                    j = i
                    while pandas.isnull(table.ix[indices[j]]['NOM p-val']):
                        j += 1
                        if j == max(indices):
                            break
                    table.loc[indices[i], 'NOM p-val'] = table['NOM p-val'][indices[j]]
                    if table.ix[indices[i]]['FDR q-val'] == 1:
                        table.loc[indices[i], 'FDR q-val'] = table['FDR q-val'][indices[j]]
                    if pandas.isnull(table.ix[indices[i]]['NES']):
                        table.loc[indices[i], 'NES'] = table['NES'][indices[j]]

        return table

    def parse_pathway_excel(self, leadingEdge=False,
                            leadingEdgeGenes=False, allGenes=False):
        self.first_pathways = None
        self.second_pathways = None

        first = glob.glob(self.directory + '/gsea_report_for*_' + self.phenotype1 + '_*xls')
        second = glob.glob(self.directory + '/gsea_report_for*_' + self.phenotype2 + '_*xls')

        assert len(first) == 1, "No file found for phenotype 1!"
        assert len(second) == 1, "No file found for phenotype 1!"

        self.first_pathways = pandas.read_table(first[0], sep='\t')
        del self.first_pathways['Unnamed: 11']

        self.first_pathways = self._check_for_nas(self.first_pathways)

        self.first_pathways = self.first_pathways.fillna('')

        self.second_pathways = pandas.read_table(second[0], sep='\t')
        del self.second_pathways['Unnamed: 11']

        self.second_pathways = self._check_for_nas(self.second_pathways)

        self.second_pathways = self.second_pathways.fillna('')

        num_in_leading_edge = []
        all_genes = []
        leading_edge_genes = []

        for pathway in self.second_pathways['NAME'].tolist():
            pathway = os.path.join(self.directory, pathway + '.xls')
            if not os.path.exists(pathway):
                n = 'NA'
                tmp_leading_edge_genes = ''
                tmp_all_genes = ''
            else:
                tmp = pandas.read_table(pathway, sep='\t')
                n = tmp[tmp['CORE ENRICHMENT'] == 'Yes'].shape[0]
                tmp_leading_edge_genes = tmp[tmp['CORE ENRICHMENT'] == 'Yes']['PROBE'].tolist()
                tmp_leading_edge_genes = ';'.join(tmp_leading_edge_genes)
                tmp_all_genes = ';'.join(tmp['PROBE'].tolist())

            num_in_leading_edge.append(n)
            leading_edge_genes.append(tmp_leading_edge_genes)
            all_genes.append(tmp_all_genes)

        if leadingEdge:
            self.second_pathways['NUM_IN_LEADING_EDGE'] = num_in_leading_edge
        if leadingEdgeGenes:
            self.second_pathways['LEADING_EDGE_GENES'] = leading_edge_genes
        if allGenes:
            self.second_pathways['ALL_GENES'] = all_genes

        num_in_leading_edge = []
        all_genes = []
        leading_edge_genes = []

        for pathway in self.first_pathways['NAME'].tolist():
            pathway = os.path.join(self.directory, pathway + '.xls')
            if not os.path.exists(pathway):
                n = 'NA'
                tmp_leading_edge_genes = ''
                tmp_all_genes = ''
            else:
                tmp = pandas.read_table(pathway, sep='\t')
                n = tmp[tmp['CORE ENRICHMENT'] == 'Yes'].shape[0]
                tmp_leading_edge_genes = tmp[tmp['CORE ENRICHMENT'] == 'Yes']['PROBE'].tolist()
                tmp_leading_edge_genes = ';'.join(tmp_leading_edge_genes)
                tmp_all_genes = ';'.join(tmp['PROBE'].tolist())

            num_in_leading_edge.append(n)
            leading_edge_genes.append(tmp_leading_edge_genes)
            all_genes.append(tmp_all_genes)

        if leadingEdge:
            self.first_pathways['NUM_IN_LEADING_EDGE'] = num_in_leading_edge
        if leadingEdgeGenes:
            self.first_pathways['LEADING_EDGE_GENES'] = leading_edge_genes
        if allGenes:
            self.first_pathways['ALL_EDGE_GENES'] = all_genes

    def write_to_excel(self, outfile, phenotype1='up', phenotype2='down'):

        assert len(phenotype1) <= 31, 'name of excel tab has to be less than 31 characters'
        assert len(phenotype2) <= 31, 'name of excel tab has to be less than 31 characters'

        workbook = xlsxwriter.Workbook(outfile)

        worksheetname = phenotype1
        sheet = workbook.add_worksheet(worksheetname)
        col = 0
        for ii in self.first_pathways.columns.tolist():
            sheet.write(0,col,ii)
            col+=1
        row = 1
        for ii in self.first_pathways.index:
            col = 0
            for jj in self.first_pathways.columns:
                sheet.write(row,col,self.first_pathways.ix[ii][jj])
                col+=1
            row+=1

        worksheetname = phenotype2
        sheet = workbook.add_worksheet(worksheetname)
        col = 0
        for ii in self.second_pathways.columns.tolist():
            sheet.write(0,col,ii)
            col+=1
        row = 1
        for ii in self.second_pathways.index:
            col = 0
            for jj in self.second_pathways.columns:
                sheet.write(row,col,self.second_pathways.ix[ii][jj])
                col+=1
            row+=1

        workbook.close()


class GSEAs():
    """
    class to store multiple GSEA results

    """

    def __init__(self,gseas=[]):
        self.gseas = gseas

    def write_all(self,outfile='',pval=1.0,FDR=1.0):

        workbook = xlsxwriter.Workbook(outfile)

        for i in range(0,len(self.gseas)):
            sheet = workbook.add_worksheet(self.gseas[i].name + ' up')

            first_pathways = self.gseas[i].first_pathways
            first_pathways = first_pathways[first_pathways['NOM p-val']<=pval]
            first_pathways = first_pathways[first_pathways['FDR q-val']<=FDR]

            col = 0
            for ii in first_pathways.columns.tolist():
                sheet.write(0,col,ii)
                col+=1
            row = 1
            for ii in first_pathways.index:
                col = 0
                for jj in first_pathways.columns:
                    sheet.write(row,col,first_pathways.ix[ii][jj])
                    col+=1
                row+=1
            #col = 0
            #row = 1
            #sheet.write(0,col,"NAME")
            #for ii in first_pathways["NAME"]:
            #    sheet.write(row,col,ii)
            #    row+=1

            sheet = workbook.add_worksheet(self.gseas[i].name + ' down')

            second_pathways = self.gseas[i].second_pathways
            second_pathways = second_pathways[second_pathways['NOM p-val']<=pval]
            second_pathways = second_pathways[second_pathways['FDR q-val']<=FDR]

            col = 0
            for ii in second_pathways.columns.tolist():
                sheet.write(0,col,ii)
                col+=1
            row = 1
            for ii in second_pathways.index:
                col = 0
                for jj in second_pathways.columns:
                    sheet.write(row,col,second_pathways.ix[ii][jj])
                    col+=1
                row+=1

            #col = 0
            #row = 1
            #sheet.write(0,col,"NAME")
            #for ii in second_pathways["NAME"]:
            #    sheet.write(row,col,ii)
            #    row+=1

        workbook.close()

    def write_unique(self,outfile='',pval=1.0,FDR=1.0):
        workbook = xlsxwriter.Workbook(outfile)
        for i in range(0,len(self.gseas)):

            #unqiue up-regulated pathways

            sheet = workbook.add_worksheet(self.gseas[i].name + ' up')

            not_i = []
            for j in range(0,len(self.gseas)):
                if i==j:
                    continue

                first_pathways = self.gseas[j].first_pathways
                first_pathways = first_pathways[first_pathways['NOM p-val']<=pval]
                first_pathways = first_pathways[first_pathways['FDR q-val']<=FDR]
                not_i += first_pathways['NAME'].tolist()

            first_pathways = self.gseas[i].first_pathways
            first_pathways = first_pathways[first_pathways['NOM p-val']<=pval]
            first_pathways = first_pathways[first_pathways['FDR q-val']<=FDR]

            unique_to_i = list(set(first_pathways['NAME'].tolist()).difference(set(not_i)))

            col = 0
            row = 1
            sheet.write(0,col,"NAME")
            for ii in unique_to_i:
                sheet.write(row,col,ii)
                row+=1

            #unqiue down-regulated pathways

            sheet = workbook.add_worksheet(self.gseas[i].name + ' down')

            not_i = []
            for j in range(0,len(self.gseas)):
                if i==j:
                    continue

                second_pathways = self.gseas[j].second_pathways
                second_pathways = second_pathways[second_pathways['NOM p-val']<=pval]
                second_pathways = second_pathways[second_pathways['FDR q-val']<=FDR]
                not_i += second_pathways["NAME"].tolist()

            second_pathways = self.gseas[i].second_pathways
            second_pathways = second_pathways[second_pathways['NOM p-val']<=pval]
            second_pathways = second_pathways[second_pathways['FDR q-val']<=FDR]

            unique_to_i = list(set(second_pathways["NAME"].tolist()).difference(set(not_i)))

            col = 0
            row = 1
            sheet.write(0,col,"NAME")
            for ii in unique_to_i:
                sheet.write(row,col,ii)
                row+=1

        workbook.close()



    def write_common(self,outfile='',pval=1.0,FDR=1.0):
        pass

    def write_compare_all(self,outfile=''):
        pass

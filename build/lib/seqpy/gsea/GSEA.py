import xlsxwriter
import pandas
import os
import numpy as np
import glob

class GSEA():
    """
    class to store and parse GSEA results

    """

    def __init__(self,directory='',name=''):
        self.directory = directory
        self.name = name

    def _check_for_nas(self,table):
        """
        Sometimes GSEA shows that significant pathways have NA for p-value
        """

        if sum(table['NOM p-val'].isnull()) > 0:
            indices = table.index.tolist()

            print 'NAs found for NOM p-val for significant pathways in:'
            print '\t\t' + self.directory

            for i in range(0,len(table.index)):
                if pandas.isnull(table.ix[indices[i]]['NOM p-val']):
                    were_there_nas = True
                    j = i
                    while pandas.isnull(table.ix[indices[j]]['NOM p-val']):
                        j+=1
                        if j == max(indices):
                            break
                    table.loc[indices[i],'NOM p-val'] = table['NOM p-val'][indices[j]]
                    if table.ix[indices[i]]['FDR q-val']==1:
                        table.loc[indices[i],'FDR q-val'] = table['FDR q-val'][indices[j]]
                    if pandas.isnull(table.ix[indices[i]]['NES']):
                        table.loc[indices[i],'NES'] = table['NES'][indices[j]]

        return table

    def parse_pathway_excel(self,reverse=False):
        self.up_pathways = None
        self.down_pathways = None

        if reverse:
            up = glob.glob(self.directory + '/gsea_report_for_na_neg_*xls')
        else:
            up = glob.glob(self.directory + '/gsea_report_for_na_pos_*xls')

        assert len(up)==1, "No up-regulated pathway files found"

        if reverse:
            down = glob.glob(self.directory + '/gsea_report_for_na_pos_*xls')
        else:
            down = glob.glob(self.directory + '/gsea_report_for_na_neg_*xls')

        assert len(down)==1, "No down-regulated pathway files found"


        self.up_pathways = pandas.read_table(up[0],sep='\t')
        del self.up_pathways['Unnamed: 11']

        self.up_pathways = self._check_for_nas(self.up_pathways)

        self.up_pathways = self.up_pathways.fillna('')

        self.down_pathways = pandas.read_table(down[0],sep='\t')
        del self.down_pathways['Unnamed: 11']

        self.down_pathways = self._check_for_nas(self.down_pathways)

        self.down_pathways = self.down_pathways.fillna('')

    def write_to_excel(self,outfile,up_tab_name='up',down_tab_name='down'):

        assert len(up_tab_name)<=31,'name of excel tab has to be less than 31 characters'
        assert len(down_tab_name)<=31,'name of excel tab has to be less than 31 characters'

        workbook = xlsxwriter.Workbook(outfile)

        worksheetname = up_tab_name
        sheet = workbook.add_worksheet(worksheetname)
        col = 0
        for ii in self.up_pathways.columns.tolist():
            sheet.write(0,col,ii)
            col+=1
        row = 1
        for ii in self.up_pathways.index:
            col = 0
            for jj in self.up_pathways.columns:
                sheet.write(row,col,self.up_pathways.ix[ii][jj])
                col+=1
            row+=1

        worksheetname = down_tab_name
        sheet = workbook.add_worksheet(worksheetname)
        col = 0
        for ii in self.down_pathways.columns.tolist():
            sheet.write(0,col,ii)
            col+=1
        row = 1
        for ii in self.down_pathways.index:
            col = 0
            for jj in self.down_pathways.columns:
                sheet.write(row,col,self.down_pathways.ix[ii][jj])
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

            up_pathways = self.gseas[i].up_pathways
            up_pathways = up_pathways[up_pathways['NOM p-val']<=pval]
            up_pathways = up_pathways[up_pathways['FDR q-val']<=FDR]

            col = 0
            for ii in up_pathways.columns.tolist():
                sheet.write(0,col,ii)
                col+=1
            row = 1
            for ii in up_pathways.index:
                col = 0
                for jj in up_pathways.columns:
                    sheet.write(row,col,up_pathways.ix[ii][jj])
                    col+=1
                row+=1
            #col = 0
            #row = 1
            #sheet.write(0,col,"NAME")
            #for ii in up_pathways["NAME"]:
            #    sheet.write(row,col,ii)
            #    row+=1

            sheet = workbook.add_worksheet(self.gseas[i].name + ' down')

            down_pathways = self.gseas[i].down_pathways
            down_pathways = down_pathways[down_pathways['NOM p-val']<=pval]
            down_pathways = down_pathways[down_pathways['FDR q-val']<=FDR]

            col = 0
            for ii in down_pathways.columns.tolist():
                sheet.write(0,col,ii)
                col+=1
            row = 1
            for ii in down_pathways.index:
                col = 0
                for jj in down_pathways.columns:
                    sheet.write(row,col,down_pathways.ix[ii][jj])
                    col+=1
                row+=1

            #col = 0
            #row = 1
            #sheet.write(0,col,"NAME")
            #for ii in down_pathways["NAME"]:
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

                up_pathways = self.gseas[j].up_pathways
                up_pathways = up_pathways[up_pathways['NOM p-val']<=pval]
                up_pathways = up_pathways[up_pathways['FDR q-val']<=FDR]
                not_i += up_pathways['NAME'].tolist()

            up_pathways = self.gseas[i].up_pathways
            up_pathways = up_pathways[up_pathways['NOM p-val']<=pval]
            up_pathways = up_pathways[up_pathways['FDR q-val']<=FDR]

            unique_to_i = list(set(up_pathways['NAME'].tolist()).difference(set(not_i)))

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

                down_pathways = self.gseas[j].down_pathways
                down_pathways = down_pathways[down_pathways['NOM p-val']<=pval]
                down_pathways = down_pathways[down_pathways['FDR q-val']<=FDR]
                not_i += down_pathways["NAME"].tolist()

            down_pathways = self.gseas[i].down_pathways
            down_pathways = down_pathways[down_pathways['NOM p-val']<=pval]
            down_pathways = down_pathways[down_pathways['FDR q-val']<=FDR]

            unique_to_i = list(set(down_pathways["NAME"].tolist()).difference(set(not_i)))

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

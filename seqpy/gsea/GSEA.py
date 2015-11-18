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

    def parse_pathway_excel(self):
        self.up_pathways = None
        self.down_pathways = None

        up = glob.glob(self.directory + '/gsea_report_for_na_pos_*xls')

        assert len(up)==1, "No up-regulated pathway files found"

        down = glob.glob(self.directory + '/gsea_report_for_na_neg_*xls')

        assert len(down)==1, "No down-regulated pathway files found"


        self.up_pathways = pandas.read_table(up[0],sep='\t')
        del self.up_pathways['Unnamed: 11']
        self.up_pathways = self.up_pathways.fillna('')

        self.down_pathways = pandas.read_table(down[0],sep='\t')
        del self.down_pathways['Unnamed: 11']
        self.down_pathways = self.down_pathways.fillna('')


class GSEAs():
    """
    class to store multiple GSEA results

    """

    def __init__(self,results=[]):
        self.results = results

    def write_all(self,outfile=''):

        workbook = xlsxwriter.Workbook(outfile)

        for comparison,data in all_data.items():
            worksheetname = comparison + '_up'
            sheet = workbook.add_worksheet(worksheetname)
            col = 0
            for ii in data['up'].columns.tolist():
                sheet.write(0,col,ii)
                col+=1
            row = 1
            for ii in data['up'].index:
                col = 0
                for jj in data['up'].columns:
                    sheet.write(row,col,data['up'].ix[ii][jj])
                    col+=1
                row+=1

            worksheetname = comparison + '_down'
            sheet = workbook.add_worksheet(worksheetname)
            col = 0
            for ii in data['down'].columns.tolist():
                sheet.write(0,col,ii)
                col+=1
            row = 1
            for ii in data['down'].index:
                col = 0
                for jj in data['down'].columns:
                    sheet.write(row,col,data['down'].ix[ii][jj])
                    col+=1
                row+=1

        workbook.close()

    def write_unique(self,outfile=''):
        pass

    def write_common(self,outfile=''):
        pass

    def write_compare_all(self,outfile=''):
        pass

import argparse

import pandas

def main(args):
    """
    Create the gct file. See the following for formatting
    http://www.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats


    filename: filename to save to. must be .gct file type
    data: pandas.DataFrame object (n by m). n being the number of
                   genes and m being the number of samples. Index must be
                   the gene names
    descriptions: vector of gene descriptions
    """

    assert args.out.find('.gct')!=-1,"Pass a .gct filename"

    data = pandas.read_csv(args.infile,index_col=0)

    samples = data.columns.tolist()
    genes = data.index.tolist()

    descriptions = ['NA']*data.shape[0]
    temp = ['NAME','DESCRIPTION'] + samples

    gct = pandas.DataFrame(columns=temp, index=data.index)
    gct['NAME'] = genes
    gct['DESCRIPTION'] = descriptions
    for i in samples:
        gct[i] = data[i]
    gct.to_csv('temp.txt',sep='\t',index=False)

    gct = open('temp.txt', 'r').read()
    fout = open(args.out, 'w')
    fout.write('#1.2\n')
    fout.write(str(data.shape[0]) + '\t' + str(len(samples)) + '\n')
    fout.write(gct)
    fout.close()


parser = argparse.ArgumentParser(description='')
parser.add_argument('--infile',type=str,help='CSV file containing expresion data',required=True)
parser.add_argument('--out',type=str,help='out file name (.gct)',required=False)
args = parser.parse_args()

main(args=args)

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


    data = pandas.read_csv(args.infile,index_col=0)
    if not args.nocls:
        print args.group1+args.group2
        data = data[args.group1+args.group2]

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
    fout = open(args.prefix+'.gct', 'w')
    fout.write('#1.2\n')
    fout.write(str(data.shape[0]) + '\t' + str(len(samples)) + '\n')
    fout.write(gct)
    fout.close()

    if not args.nocls:
        fout = open(args.prefix + '.cls','w')
        fout.write(str(len(args.group1) + len(args.group2)) + ' 2 1\n')
        fout.write('# ' + args.phenotypes[0] + ' ' + args.phenotypes[1] + '\n')
        fout.write(' '.join(['0']*len(args.group1)) + ' ' + ' '.join(['1']*len(args.group2)) + '\n')
        fout.close()


parser = argparse.ArgumentParser(description='Creates a gct and cls file from a expression file')
parser.add_argument('--infile',type=str,help='CSV file containing expresion data',required=True)
parser.add_argument(
    '--nocls',
    action='store_true',
    help='Do no produce a CLS file',
    required=False
)
parser.add_argument(
    '--group1',
    type=str,
    required=False,
    nargs='+',
    help='Space-delimited list of sample names'
)
parser.add_argument(
    '--group2',
    type=str,
    required=False,
    nargs='+',
    help='Space-delimited list of sample names'
)
parser.add_argument(
    '--phenotypes',
    type=str,
    required=False,
    nargs='+',
    help='Space separated name for the groups/phenotype (e.g. WT,MUT), same order and group1 and group2'
)
parser.add_argument('--prefix',type=str,help='out file prefix',required=False)
args = parser.parse_args()

main(args=args)

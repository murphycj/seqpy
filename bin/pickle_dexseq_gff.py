import sys
import re
import argparse
import pickle

def main(args):

    fin = open(args.gff,'r')
    current_gene = ''

    gene_data = dict()
    exon_intervals = dict()

    for line in fin.readlines():
        line = line.rstrip().split('\t')

        if line[2]=='aggregate_gene':
            tmp = line[8]
            tmp = tmp.replace("gene_id \"","")
            tmp = tmp.replace("\"","")
            current_gene = tmp

            gene_data[current_gene] = {}
        else:
            exon_number = re.findall("exonic_part_number \"(.*)\";", line[8])
            assert len(exon_number)==1, "error parsing exon number: " + str(exon_number)
            exon_number = exon_number[0]

            transcripts = re.findall("transcripts \"(.*)\"; exonic_part_number", line[8])
            assert len(transcripts)==1, "error parsing exon number: " + str(transcripts)
            transcripts = transcripts[0].split('+')

            exon_intervals[current_gene + ':' + exon_number] = [line[0], line[3], line[4]]

            for transcript in transcripts:
                if transcript not in gene_data[current_gene]:
                    gene_data[current_gene][transcript] = [exon_number]
                else:
                    gene_data[current_gene][transcript].append(exon_number)

    fin.close()

    data = {
        'gene_data': gene_data,
        'exon_intervals': exon_intervals
    }

    fout = open(args.gff.replace('.gff','.pickle'),'wb')
    pickle.dump(data,fout)
    fout.close()

parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--gff',
    type=str,
    required=True,
    help='The DEXseq gff file.'
)
args = parser.parse_args()

main(args=args)

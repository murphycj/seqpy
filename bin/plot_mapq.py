import argparse
import pysam
import plotnine as pn
import pandas as pd

def main(args):

    bamin = bamin = pysam.AlignmentFile(open(args.i,'rb'))
    mapqs = []

    for record in bamin:
        if not record.is_unmapped:
            mapqs.append(record.mapq)

    mapqs = pd.DataFrame(mapqs)
    mapqs.columns = ['MAPQ']

    p = pn.ggplot(mapqs,pn.aes('MAPQ')) + pn.geom_histogram()
    p.save(args.o)

parser = argparse.ArgumentParser(description='')
parser.add_argument('--i',type=str,help='In file name',required=True)
parser.add_argument('--o',type=str,help='Output file name',required=True)
args = parser.parse_args()

main(args=args)

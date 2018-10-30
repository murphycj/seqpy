import argparse
import vcf
import plotnine as pn
import pandas as pd

def main(args):

    vin = vcf.Reader(open(args.i,'r'))
    data = {
        'vaf':[],
        'DP':[]
    }

    for variant in vin:
        ad = variant.samples[0].data.AD

        if len(ad)!=2:
            continue

        if ad[0]==0 or ad[1]==0:
            continue

        data['vaf'].append(float(ad[1])/sum(ad))
        data['DP'].append(variant.INFO['DP'])

    data = pd.DataFrame(data)
    data = data[data['DP']>9]

    p = pn.ggplot(data,pn.aes('vaf')) + pn.geom_histogram(binwidth=0.01) + pn.xlim(0,1)
    p.save(args.o)

    p = pn.ggplot(data,pn.aes(x='vaf',y='DP')) + pn.geom_point() + pn.xlim(0,1)
    p.save(args.o2)


parser = argparse.ArgumentParser(description='')
parser.add_argument('--i',type=str,help='In file name',required=True)
parser.add_argument('--o',type=str,help='Output file name',required=True)
parser.add_argument('--o2',type=str,help='Output file name',required=True)
args = parser.parse_args()

main(args=args)

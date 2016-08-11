import pandas
import argparse

def main(args):
    data = zip(args.files,args.samples)
    results = {}
    samples = []
    locations = []

    for d in data:
        snpden = pandas.read_table(d[0], sep='\t', index_col=None,low_memory=False)
        temp = zip(snpden['CHROM'],snpden['BIN_START'])
        temp2 = pandas.DataFrame(snpden['VARIANTS/KB'])
        temp2.index = temp
        temp2.columns = [d[1]]
        results[d[1]] = temp2
        samples.append(d[1])
        locations += temp

    locations = list(set(locations))

    all_data = pandas.DataFrame(index = locations, columns = ['CHROM','BIN_START'] + samples)
    all_data['CHROM'] = map(lambda x: x[0],locations)
    all_data['BIN_START'] = map(lambda x: x[1],locations)

    for s, d in results.items():
        all_data.loc[:,s] = d

    all_data=all_data.fillna(0)

    temp = all_data[samples]>=args.min_den
    all_data = all_data[temp.sum(axis=1)>0]
    all_data = all_data.sort_index()
    all_data.to_csv(args.out,index=False)


parser = argparse.ArgumentParser(description='Aggregates the snpden values from cufflinks output')
parser.add_argument(
    '--files',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list.'
)
parser.add_argument(
    '--samples',
    type=str,
    required=True,
    nargs='+',
    help='Space-delimited list of sample names'
)
parser.add_argument(
    '--min_den',
    type=int,
    help='Filter for only windows where at least one samples has min_den (default 5)',
    default=5
)
parser.add_argument(
    '--out',
    type=str,
    required=True,
    help='Output file name'
)
args = parser.parse_args()

main(args=args)

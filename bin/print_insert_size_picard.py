import argparse
from seqpy.parsers import CollectInsertSizeMetrics


def main(args):
    insert_data = CollectInsertSizeMetrics(args.file)

    if args.round:
        print int(insert_data.metrics[args.type])
    else:
        print insert_data.metrics[args.type]


parser = argparse.ArgumentParser(
    description='Prints the median insert size from Picardtools ' +
                'CollectInsertSizeMetrics'
)
parser.add_argument(
    '--file',
    type=str,
    required=True,
    help='Picard output file containing the insert size metrics.'
)
parser.add_argument(
    '--type',
    type=str,
    required=False,
    default='MEAN_INSERT_SIZE',
    help='(optional) What metric to print (default MEAN_INSERT_SIZE)'
)
parser.add_argument(
    '--round',
    required=False,
    action='store_true'
    help='Round to nearest integer'
)
args = parser.parse_args()

main(args=args)

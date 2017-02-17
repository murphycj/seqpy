import argparse
from seqpy.parsers import CollectInsertSizeMetrics


def main(args):
    insert_data = CollectInsertSizeMetrics(args.file)

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
    required=Flase,
    default='MEAN_INSERT_SIZE',
    help='(optional) What metric to print (default MEAN_INSERT_SIZE)'
)
args = parser.parse_args()

main(args=args)

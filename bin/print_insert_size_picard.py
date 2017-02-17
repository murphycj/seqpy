import argparse
from seqpy.parsers import CollectInsertSizeMetrics


def main(args):
    insert_data = CollectInsertSizeMetrics(args.file)

    print insert_data.metrics['MEAN_INSERT_SIZE']


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
args = parser.parse_args()

main(args=args)

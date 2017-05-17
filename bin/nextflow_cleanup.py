import os
import argparse
import glob

def main(args):
    dirs_to_delete = []
    for f in glob.glob(args.work + '/*'):
        for ff in glob.glob(f + '/*'):
            exit_file = os.path.join(ff,'.exitcode')
            if os.path.exists(exit_file):
                exitcode = open(exit_file,'r').read()
                if int(exitcode) not in args.exit_code:
                    dirs_to_delete.append(ff)
    for i in dirs_to_delete:
        cmd = 'rm -rf ' + os.path.abspath(i)
        print cmd
        os.system(cmd)

parser = argparse.ArgumentParser(description='Remove nextflow work directories that failed')
parser.add_argument(
    '-w',
    '--work',
    type=str,
    required=True,
    help='Nextflow work directory'
)
parser.add_argument(
    '-e',
    '--exit_code',
    type=int,
    required=False,
    nargs='+',
    default=[0],
    help='Space-delimited list of exit codes to not remove.'
)
args = parser.parse_args()

main(args=args)

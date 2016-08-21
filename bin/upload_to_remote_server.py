"""
Copy a directory to a remote server via lftp
"""

import argparse
import os, glob

def main(args):

    #create root directory on server

    if not args.remoteroot==".":
        cmd = + args.lftp + ' ' + \
            '-c \"open -u ' + args.user + ',' + args.password + ' ' + \
            args.remote + '; ' + \
            'mkdir -f ~/' + args.remoteroot + ';\"'
        print cmd

    #copy directory structure

    for d in os.walk(args.dir):
        cmd = + args.lftp + ' ' + \
            '-c \"open -u ' + args.user + ',' + args.password + ' ' + \
            args.remote + '; ' + \
            'mkdir -f ~/' + args.remoteroot + '/' + \
            d[0] + ';\"'
        print cmd

    #upload files

    for d in os.walk(args.dir):
        for f in d[2]:
            cmd = + args.lftp + ' ' + \
                '-c \"open -u ' + args.user + ',' + args.password + ' ' + \
                args.remote + '; ' + \
                'put -O ~/' + args.remoteroot + '/' + d[0] + '/' + f + ';\"'
            print cmd


parser = argparse.ArgumentParser(description='')
parser.add_argument(
    '--dir',
    type=str,
    required=True,
    help='Directory to copy to remote server'
)
parser.add_argument(
    '--remoteroot',
    type=str,
    required=True,
    help='Root directory to store results on server',
)
parser.add_argument(
    '--remote',
    type=str,
    default='sftp://serapeum1.qib.pbtech',
    help='Remote server address (default sftp://serapeum1.qib.pbtech)'
)
parser.add_argument(
    '--lftp',
    type=str,
    require=True,
    help='Path to lftp'
)
parser.add_argument(
    '--user',
    type=str,
    required=True
)
parser.add_argument(
    '--password',
    type=str,
    required=True
)
args = parser.parse_args()

main(args=args)

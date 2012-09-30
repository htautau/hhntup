#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-u', '--user')
parser.add_argument('-d', '--dir')
parser.add_argument('-v', '--verbose', action='store_true', default=False)
parser.add_argument('files', nargs='*')
args = parser.parse_args()

from glob import glob
import os
import subprocess

AFS_WORKSPACE = '/afs/cern.ch/work/%s/%s/public' % (args.user[0], args.user)
HOST = '%s@lxplus5.cern.ch' % args.user

if args.dir:
    AFS_WORKSPACE = os.path.join(AFS_WORKSPACE, args.dir)

def cmd(args, verbose=args.verbose):

    if verbose:
        print ' '.join(args)
    subprocess.call(args)

cmd(['scp'] + args.files + [':'.join([HOST, AFS_WORKSPACE])])

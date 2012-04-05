#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-u', '--user')
parser.add_argument('dir')
args = parser.parse_args()

from glob import glob
import os
import subprocess

AFS_WORKSPACE = '/afs/cern.ch/work/%s/%s/public' % (args.user[0], args.user)
HOST = '%s@lxplus5.cern.ch' % args.user

subprocess.call(['scp'] + glob(os.path.join(args.dir, '*.root')) + [':'.join([HOST, AFS_WORKSPACE])])

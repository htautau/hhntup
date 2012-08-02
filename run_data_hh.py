#!/usr/bin/env python

import os
import cluster

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--nproc', type=int, default=3)
parser.add_argument('--nsplit', type=int, default=30)
parser.add_argument('--queue', default='short')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--dry', action='store_true', default=False)
args = parser.parse_args()

setup = cluster.get_setup('setup.noel.sfu.txt')

CWD = os.getcwd()
CMD = ("%s && ./run --output-path ntuples/hadhad "
       "-s HHProcessor.py -n %d --db datasets_hh "
       "--nice %d --split %d:%%d data") % (
               setup, args.nproc, args.nice, args.nsplit)

for i in xrange(args.nsplit):
    cmd = "cd %s && %s" % (CWD, CMD % (i + 1))
    cluster.qsub(
        cmd,
        ncpus=args.nproc,
        name='HHProcessor.data_%d' % (i + 1),
        stderr_path='ntuples/hadhad',
        stdout_path='ntuples/hadhad',
        queue=args.queue,
        dry_run=args.dry)

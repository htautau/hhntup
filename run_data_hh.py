#!/usr/bin/env python

import os
import cluster

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--year', type=int, choices=(11, 12), default=11)
parser.add_argument('--student', default='HHProcessor.py')
parser.add_argument('--nproc', type=int, default=3)
parser.add_argument('--nsplit', type=int, default=30)
parser.add_argument('--queue', default='short')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--dry', action='store_true', default=False)
parser.add_argument('--output-path',
                    default='/cluster/data11/endw/ntuples/hadhad_running')
parser.add_argument('--db', default='datasets_hh')
parser.add_argument('splits', nargs='*', type=int)
args = parser.parse_args()

setup = cluster.get_setup('setup.noel.sfu.txt')

student_name = os.path.splitext(args.student)[0]

output_path = os.path.join(args.output_path, student_name)

CWD = os.getcwd()
CMD = ("%s && ./run --output-path %s "
       "-s %s -n %d --db %s "
       "--nice %d --split %d:%%d data%d-JetTauEtmiss") % (
               setup, output_path, args.student, args.nproc,
               args.db, args.nice, args.nsplit, args.year)

for i in xrange(args.nsplit):
    if args.splits and (i + 1) not in args.splits:
        continue
    cmd = "cd %s && %s" % (CWD, CMD % (i + 1))
    output = '%s.data%d-JetTauEtmiss_%d.root' % (student_name, args.year, i + 1)
    if os.path.exists(os.path.join(output_path, output)):
        print "%s already exists. please delete it and resubmit" % output
        continue
    cluster.qsub(
        cmd,
        ncpus=args.nproc,
        name='%s.data%d_%d' % (student_name, args.year, i + 1),
        stderr_path=output_path,
        stdout_path=output_path,
        queue=args.queue,
        dry_run=args.dry)

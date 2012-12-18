#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--year', type=int, choices=(11, 12), default=11)
parser.add_argument('--student', default='HHProcessor.py')
parser.add_argument('--systematics', default=None)
parser.add_argument('--nproc', type=int, default=5)
parser.add_argument('--nsplit', type=int, default=30)
parser.add_argument('--queue', default='short')
parser.add_argument('--output-path', default='ntuples/hadhad_running')
parser.add_argument('--db', default='datasets_hh')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--nominal-only', action='store_true', default=False)
parser.add_argument('--systematics-only', action='store_true', default=False)
parser.add_argument('--dry', action='store_true', default=False)
parser.add_argument('--use-ssh', dest='use_qsub', action='store_false', default=True)
parser.add_argument('--warnings-as-errors', action='store_true', default=False)
parser.add_argument('splits', nargs='*', type=int)

args = parser.parse_args()

import sys
import os
import cluster
from higgstautau import samples

hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

student_name = os.path.splitext(args.student)[0]
output_path = os.path.join(args.output_path, student_name)

CWD = os.getcwd()
CMD = ("%s && ./run --output-path %s "
       "-s %s -n %d --db %s "
       "--nice %d --split %d:%%d %%s") % (
               setup, output_path,
               args.student,
               args.nproc,
               args.db,
               args.nice,
               args.nsplit)

def run_sample(sample):

    for i in xrange(args.nsplit):
        if args.splits and (i + 1) not in args.splits:
            continue
        cmd = "cd %s && %s" % (CWD, CMD % (i + 1, sample))
        job_name = '%s.%s_%d' % (student_name, sample, i + 1)
        output = job_name + '.root'
        if os.path.exists(os.path.join(output_path, output)):
            print "%s already exists. please delete it and resubmit" % output
            continue
        cluster.qsub(
            cmd,
            ncpus=args.nproc,
            name=job_name,
            stderr_path=output_path,
            stdout_path=output_path,
            queue=args.queue,
            dry_run=args.dry)

if not args.systematics_only:
    # nominal values
    datasets = samples.samples('hadhad', args.year, '*embed*')
    for sample in datasets:
        run_sample(sample)

if not args.nominal_only:
    if args.systematics is not None:
        args.systematics = [
                tuple(s.upper().split('+')) for s in
                args.systematics.split(',')]
    # systematics
    for datasets, systematics in samples.iter_samples('hadhad', args.year,
            '*embed*', systematics=True):
        print datasets, systematics

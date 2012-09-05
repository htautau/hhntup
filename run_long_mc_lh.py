#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--student', default='LHProcessor.py')
parser.add_argument('--systematics', default=None)
parser.add_argument('--nproc', type=int, default=12)
parser.add_argument('--queue', default='medium')
parser.add_argument('--output-path', default='ntuples/lephad')
parser.add_argument('--db', default='datasets_lh')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--nominal-only', action='store_true', default=False)
parser.add_argument('--systematics-only', action='store_true', default=False)
parser.add_argument('--dry', action='store_true', default=False)
parser.add_argument('--use-ssh', dest='use_qsub', action='store_false', default=True)
parser.add_argument('samples', nargs='?', default=None)

args = parser.parse_args()

import sys
import cluster
from higgstautau import samples

hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.michel.sfu.txt')

datasets = [
    'T1_McAtNlo_Jimmy.mc11c',
    'AlpgenJimmyZeeNp0_pt20.mc11c',
    'AlpgenJimmyZtautauNp0_pt20.mc11c',
    ]

if not args.systematics_only:
    # nominal values
    cluster.run(args.student,
                db=args.db,
                datasets=datasets,
                hosts=hosts,
                nproc=args.nproc,
                nice=args.nice,
                setup=setup,
                output_path=args.output_path,
                use_qsub=args.use_qsub,
                qsub_queue=args.queue,
                dry_run=args.dry)

if not args.nominal_only:
    if args.systematics is not None:
        args.systematics = [
                set(s.upper().split('+')) for s in
                args.systematics.split(',')]
    # systematics
    cluster.run_systematics('LEPHAD',
                args.student,
                db=args.db,
                systematics=args.systematics,
                datasets=datasets,
                hosts=hosts,
                nproc=args.nproc,
                nice=args.nice,
                setup=setup,
                output_path=args.output_path,
                use_qsub=args.use_qsub,
                qsub_queue=args.queue,
                dry_run=args.dry)

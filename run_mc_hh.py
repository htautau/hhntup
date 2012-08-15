#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--systematics', default=None)
parser.add_argument('--nproc', type=int, default=5)
parser.add_argument('--queue', default='short')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--nominal-only', action='store_true', default=False)
parser.add_argument('--systematics-only', action='store_true', default=False)
parser.add_argument('--dry', action='store_true', default=False)
parser.add_argument('samples', nargs='?', default=None)

args = parser.parse_args()

import sys
import cluster
from higgstautau import samples

hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

datasets = samples.samples('hadhad', args.samples)
if not datasets:
    sys.exit('No datasets selected!')

if not args.systematics_only:
    # nominal values
    cluster.run('HHProcessor.py',
                db='datasets_hh',
                datasets=datasets,
                hosts=hosts,
                nproc=args.nproc,
                nice=args.nice,
                setup=setup,
                output_path='ntuples/hadhad',
                use_qsub=True,
                qsub_queue=args.queue,
                dry_run=args.dry)

if not args.nominal_only:
    if args.systematics is not None:
        args.systematics = [
                set(s.upper().split('+')) for s in
                args.systematics.split(',')]
    # systematics
    cluster.run_systematics('HADHAD',
                'HHProcessor.py',
                db='datasets_hh',
                systematics=args.systematics,
                datasets=datasets,
                hosts=hosts,
                nproc=args.nproc,
                nice=args.nice,
                setup=setup,
                output_path='ntuples/hadhad',
                use_qsub=True,
                qsub_queue=args.queue,
                dry_run=args.dry)

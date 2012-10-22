#!/usr/bin/env python

import os
import cluster

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--nproc', type=int, default=6)
parser.add_argument('--nsplit', type=int, default=30)
parser.add_argument('--queue', default='short')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--dry', action='store_true', default=False)
args = parser.parse_args()

setup = cluster.get_setup('setup.michel.sfu.txt')

datasets = [#'data-Egamma',
            #'data-JetTauEtmiss',
            #'data-Muons',
            'embed-lh-isol-mfsim'
            ]

CWD = os.getcwd()

for i in xrange(args.nsplit):
    for dataset in datasets:

        CMD = ("%s && ./run --output-path ntuples/lephadComposite/%s "
               "-s LHProcessor.py -n %d --db datasets_lh "
        "--nice %d --split %d:%%d %s") % (
            setup, dataset, args.nproc, args.nice, args.nsplit, dataset)
        
        cmd = "cd %s && %s" % (CWD, CMD % (i + 1))
        cluster.qsub(
            cmd,
            ncpus=args.nproc,
            name= ('LHProcessor.data-%s_%d') % (dataset, (i + 1)),
            stderr_path='ntuples/lephadComposite/' + dataset,
            stdout_path='ntuples/lephadComposite/' + dataset,
            queue=args.queue,
            dry_run=args.dry)

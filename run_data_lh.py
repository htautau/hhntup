#!/usr/bin/env python

import os
import cluster

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--nproc', type=int, default=5)
parser.add_argument('--nsplit', type=int, default=20)
parser.add_argument('--queue', default='short')
parser.add_argument('--nice', type=int, default=10)
parser.add_argument('--dry', action='store_true', default=False)
args = parser.parse_args()

setup = cluster.get_setup('setup.michel.sfu.txt')

datasets2011 = [
    'data11-Egamma',
    'data11-JetTauEtmiss',
    'data11-Muons',
    'embed11-lh-isol-mfsim'
    ]

datasets2012 = [
    #'data12-Egamma',
    'data12-JetTauEtmiss',
    'data12-PhysCont',
    #'embed12-lh-mfsim',
    #'embed12-lh-mfsup',
    #'embed12-lh-mfsdn'
    ]

datasets = datasets2012

CWD = os.getcwd()

for i in xrange(args.nsplit):
    for dataset in datasets:

        CMD = ("%s && ./run --output-path /cluster/data05/michel/Ntuples/lephad/2012/%s "
               "-s LHProcessorCN.py -n %d --db datasets_lh "
        "--nice %d --split %d:%%d %s") % (
            setup, dataset, args.nproc, args.nice, args.nsplit, dataset)
        
        cmd = "cd %s && %s" % (CWD, CMD % (i + 1))
        cluster.qsub(
            cmd,
            ncpus=args.nproc,
            name= ('LHProcessorCN.data-%s_%d') % (dataset, (i + 1)),
            stderr_path='/cluster/data05/michel/Ntuples/lephad/2012/' + dataset,
            stdout_path='/cluster/data05/michel/Ntuples/lephad/2012/' + dataset,
            queue=args.queue,
            dry_run=args.dry)

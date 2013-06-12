#!/usr/bin/env python

import os
import cluster

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('--nproc', type=int, default=5)
parser.add_argument('--nsplit', type=int, default=1)
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
    #'embed12-LH-IM',
    'data12-PhysCont-periodJ.JetTauEtmiss',
    #'data12-PhysCont-periodH.Egamma',
    #'data12-PhysCont-periodJ.Muons',
    #'data12-PhysCont-periodH.Muons',
    #'data12-PhysCont-periodA.Muons',
    #'data12-PhysCont-periodB.Muons',
    #'data12-PhysCont-periodD.Muons',
    #'data12-PhysCont-periodA.Egamma',
    #'data12-PhysCont-periodE.JetTauEtmiss',
    #'data12-PhysCont-periodI.JetTauEtmiss',
    #'data12-PhysCont-periodH.JetTauEtmiss',
    #'data12-PhysCont-periodG.Egamma',
    #'data12-PhysCont-periodC.Egamma',
    #'data12-PhysCont-periodL.Egamma',
    #'data12-PhysCont-periodB.JetTauEtmiss',
    #'data12-PhysCont-periodG.JetTauEtmiss',
    #'data12-PhysCont-periodE.Muons',
    #'data12-PhysCont-periodD.JetTauEtmiss',
    #'data12-PhysCont-periodC.JetTauEtmiss',
    #'data12-PhysCont-periodC.Muons',
    #'data12-PhysCont-periodL.JetTauEtmiss',
    #'data12-PhysCont-periodG.Muons',
    #'data12-PhysCont-periodL.Muons',
    #'data12-PhysCont-periodE.Egamma',
    #'data12-PhysCont-periodD.Egamma',
    #'data12-PhysCont-periodA.JetTauEtmiss',
    #'data12-PhysCont-periodB.Egamma',
    #'data12-PhysCont-periodI.Muons',
    #'data12-PhysCont-periodI.Egamma',
    #'data12-PhysCont-periodJ.Egamma',
    ]

datasets = datasets2012

CWD = os.getcwd()

for i in xrange(args.nsplit):
    for dataset in datasets:

        CMD = ("%s && ./run --output-path /cluster/data05/michel/Ntuples/lephad/2012-p1443-4 "
               "-s LHProcessorCN.py -n %d --db datasets_lh "
        "--nice %d --split %d:%%d %s") % (
            setup, args.nproc, args.nice, args.nsplit, dataset)
        
        cmd = "cd %s && %s" % (CWD, CMD % (i + 1))
        cluster.qsub(
            cmd,
            ncpus=args.nproc,
            name= ('LHProcessorCN.data-%s_%d') % (dataset, (i + 1)),
            stderr_path='/cluster/data05/michel/Ntuples/lephad/2012-p1443-4',
            stdout_path='/cluster/data05/michel/Ntuples/lephad/2012-p1443-4',
            queue=args.queue,
            dry_run=args.dry)

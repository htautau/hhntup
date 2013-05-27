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
    # 'data12-v03',
    # 'data12-v03-periodA.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodA.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodA.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodB.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodB.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodB.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodC.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodC.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodC.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodD.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodD.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodD.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodE.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodE.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodE.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodG.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodG.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodG.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodH.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodH.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodH.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodI.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodI.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodI.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodJ.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodJ.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodJ.Muons.TunaCont.2013-March-29',
    # 'data12-v03-periodL.Egamma.TunaCont.2013-March-29',
    # 'data12-v03-periodL.JetTauEtmiss.TunaCont.2013-March-29',
    # 'data12-v03-periodL.Muons.TunaCont.2013-March-29',
    'embed12-LH-DN',
    'embed12-LH-IM',
    'embed12-LH-UP',
    ]

datasets = datasets2012

CWD = os.getcwd()

for i in xrange(args.nsplit):
    for dataset in datasets:

        CMD = ("%s && ./run --output-path /cluster/data05/michel/Ntuples/lephad/2012-p1443-2 "
               "-s LHProcessorCN.py -n %d --db datasets_lh "
        "--nice %d --split %d:%%d %s") % (
            setup, args.nproc, args.nice, args.nsplit, dataset)
        
        cmd = "cd %s && %s" % (CWD, CMD % (i + 1))
        cluster.qsub(
            cmd,
            ncpus=args.nproc,
            name= ('LHProcessorCN.data-%s_%d') % (dataset, (i + 1)),
            stderr_path='/cluster/data05/michel/Ntuples/lephad/2012-p1443-2',
            stdout_path='/cluster/data05/michel/Ntuples/lephad/2012-p1443-2',
            queue=args.queue,
            dry_run=args.dry)

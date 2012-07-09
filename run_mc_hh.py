#!/usr/bin/env python

import sys
import cluster
from higgstautau import samples

hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

datasets = samples.samples('hadhad')

# nominal values
cluster.run('HHProcessor.py',
            db='datasets_hh',
            datasets=datasets,
            hosts=hosts,
            nproc=5,
            nice=10,
            setup=setup,
            output_path='ntuples/hadhad',
            use_qsub=True,
            qsub_queue='short',
            dry_run='dry' in sys.argv)

# systematics
cluster.run_systematics('HADHAD',
            'HHProcessor.py',
            db='datasets_hh',
            datasets=datasets,
            hosts=hosts,
            nproc=5,
            nice=10,
            setup=setup,
            output_path='ntuples/hadhad',
            use_qsub=True,
            qsub_queue='short',
            dry_run='dry' in sys.argv)

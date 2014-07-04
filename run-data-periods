#!/usr/bin/env python

import cluster


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

datasets = [
    'data-JetTauEtmiss-B',
    'data-JetTauEtmiss-D',
    'data-JetTauEtmiss-E',
    'data-JetTauEtmiss-F',
    'data-JetTauEtmiss-G',
    'data-JetTauEtmiss-H',
    'data-JetTauEtmiss-I',
    'data-JetTauEtmiss-J',
    'data-JetTauEtmiss-K',
    'data-JetTauEtmiss-L',
    'data-JetTauEtmiss-M',
]

cluster.run('HHProcessor.py',
            db='datasets_hh',
            datasets=datasets,
            hosts=hosts,
            nproc=10,
            nice=10,
            setup=setup,
            use_qsub=True)

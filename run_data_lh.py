#!/usr/bin/env python

import cluster


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.michel.sfu.txt')

datasets = [
'data'
]

cluster.run('muLHProcessor.py',
            datasets=datasets,
            hosts=hosts,
            nproc=8,
            nice=10,
            setup=setup)

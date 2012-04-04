#!/usr/bin/env python

import cluster


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.txt')

datasets = [
    'data-B',
    'data-D',
    'data-E',
    'data-F',
    'data-G',
    'data-H',
    'data-I',
    'data-J',
    'data-K',
    'data-L',
    'data-M',
]

cluster.run('HHProcessor.py',
            datasets=datasets,
            hosts=hosts,
            nproc=10,
            nice=10,
            setup=setup)

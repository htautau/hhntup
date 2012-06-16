#!/usr/bin/env python

import sys
import cluster


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

cluster.run('HHProcessor.py',
            db='datasets_hh',
            datasets=sys.argv[1:],
            hosts=hosts,
            nproc=10,
            nice=10,
            setup=setup)

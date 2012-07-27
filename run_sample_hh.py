#!/usr/bin/env python

import sys
import cluster


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

cluster.run('HHProcessor.py',
            db='datasets_hh',
            datasets=[sys.argv[1]],
            hosts=hosts,
            nproc=1,
            nice=10,
            setup=setup,
            output_path='.',
            student_args=sys.argv[2:],
            use_qsub=False)

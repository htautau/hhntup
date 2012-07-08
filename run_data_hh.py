#!/usr/bin/env python

import os
import cluster


setup = cluster.get_setup('setup.noel.sfu.txt')

CWD = os.getcwd()
NPROC = 5
NSPLIT = 30
CMD = "%s && ./run -s HHProcessor.py -n %d --db datasets_hh --nice 10 --split %d:%%d data" % (setup, NPROC, NSPLIT)

for i in xrange(1, NSPLIT + 1):
    cmd = "cd %s && %s" % (CWD, CMD % i)
    cluster.qsub(cmd, ncpus=NPROC)

#!/usr/bin/env python

import os
import subprocess
import socket
import shlex
import multiprocessing as mp
from itertools import cycle


HOSTNAME = socket.gethostname()
CWD = os.getcwd()

hosts = ['lhc%02d' % i for i in xrange(1, 11)]
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

NPROC = 10
CMD = "./run -s HHProcessor.py -n %d --nice 10 %%s" % NPROC

proc_cmds = []

setup = 'source /cluster/data10/endw/bashrc.sfu/bashrc'

for host in cycle(hosts):
    if len(datasets) == 0:
        break
    ds = datasets.pop(0)
    cmd = CMD % ds
    print host
    if not HOSTNAME.startswith(host):
        cmd = "ssh %s '%s && cd %s && %s'" % (host, setup, CWD, cmd)
    print cmd
    proc_cmds.append(cmd)


def run(cmd):

    subprocess.call(cmd, shell=True)


procs = []
for cmd in proc_cmds:
    proc = mp.Process(target=run, args=(cmd,))
    proc.start()
    procs.append(proc)
for proc in procs:
    proc.join()

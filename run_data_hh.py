#!/usr/bin/env python

import os
import subprocess
import socket
import shlex
import multiprocessing as mp


HOSTNAME = socket.gethostname()
CWD = os.getcwd()

hosts = ['lhc%02d' % i for i in xrange(1, 11)]

NPROC = 10
CMD = "./run -s HHProcessor.py -n %d --nice 10 --split %d:%%d data" % (NPROC, len(hosts))

proc_cmds = []

setup = 'source /cluster/data10/endw/bashrc.sfu/bashrc'

for i, host in enumerate(hosts):
    cmd = CMD % (i + 1)
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

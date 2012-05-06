#!/usr/bin/env python

import cluster
import socket
import os
import multiprocessing as mp
import subprocess


hosts = cluster.get_hosts('hosts.sfu.txt')
setup = cluster.get_setup('setup.noel.sfu.txt')

HOSTNAME = socket.gethostname()
CWD = os.getcwd()

NPROC = 10
CMD = "%s && ./run -s HHProcessor.py -n %d --db datasets_hh --nice 10 --split %d:%%d data" % (setup, NPROC, len(hosts))

proc_cmds = []

for i, host in enumerate(hosts):
    cmd = CMD % (i + 1)
    cmd = "ssh %s 'cd %s && %s'" % (host, CWD, cmd)
    print "%s: %s" % (host, cmd)
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

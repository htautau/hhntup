#!/usr/bin/env python

import os
import subprocess
import socket
import shlex
import multiprocessing as mp

tasks = {
    'lhc07': ["PowHegPythia_VBFH*_tautauhh.mc11c"],
    'lhc10': ["PowHegPythia_ggH*_tautauhh.mc11c"],
    'lhc08': ["AlpgenJimmyZ*Np[0-5]_pt20.mc11c"],
    'lhc09': ["AlpgenJimmyW*Np[0-5]_pt20.mc11c"],
}

HOSTNAME = socket.gethostname()
CWD = os.getcwd()
NPROC = 4
NICE = 10
CMD = "./run -s HHProcessor.py --nproc %d --nice %d " % (NPROC, NICE)

proc_cmds = []

setup = 'source /cluster/data10/endw/bashrc.sfu/bashrc'

for host, samples in tasks.items():
    cmd = CMD + " ".join(['"%s"' % s for s in samples])
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

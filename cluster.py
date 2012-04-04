import socket
import os
import subprocess
import shlex
import multiprocessing as mp
from itertools import cycle
import shlex


HOSTNAME = socket.gethostname()


def get_load(host):

    cmd = 'python -c "import os; print os.getloadavg()[1]"'
    if not HOSTNAME.startswith(host):
        cmd = "ssh %s '%s'" % (host, cmd)
    load = float(subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].strip())
    return load


def get_hosts(filename):

    hosts = None
    with open(filename) as f:
        hosts = [line.strip() for line in f.readlines()]
    return hosts


def get_setup(filename):

    with open(filename) as f:
        return ' && '.join([line.strip() for line
                            in f.readlines()])


def run(student, datasets, hosts,
        nproc=1,
        nice=0,
        setup=None):

    CMD = "./run -s %s -n %d --nice %d %%s" % (student, nproc, nice)
    CWD = os.getcwd()

    proc_cmds = []
    for host in cycle(hosts):
        if len(datasets) == 0:
            break
        ds = datasets.pop(0)
        cmd = CMD % ds
        print host
        if not HOSTNAME.startswith(host):
            if setup is not None:
                cmd = "ssh %s '%s && cd %s && %s'" % (host, setup, CWD, cmd)
            else:
                cmd = "ssh %s 'cd %s && %s'" % (host, CWD, cmd)
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

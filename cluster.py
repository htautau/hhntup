import socket
import os
import subprocess
import multiprocessing as mp


HOSTNAME = socket.gethostname()


class Host(object):

    def __init__(self, name):

        self.name = name
        self.njobs = 0

    @property
    def load(self):

        return get_load(self.name)

    def load_metric(self):

        return self.load * (self.njobs + 1) + self.njobs

    def __cmp__(self, other):

        return cmp(self.load_metric(),
                   other.load_metric())

    def __str__(self):

        return "%s(%.3f:%d)" % (self.name, self.load, self.njobs)

    def __repr__(self):

        return str(self)


def get_load(host):

    cmd = 'python -c "import os; print os.getloadavg()[0]"'
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


def run_helper(cmd):

        subprocess.call(cmd, shell=True)


def run(student,
        db,
        datasets,
        hosts,
        nproc=1,
        nice=0,
        setup=None):

    CMD = "./run -s %s -n %d --db %s --nice %d %%s" % (student, nproc, db, nice)
    if setup is not None:
        CMD = "%s && %s" % (setup, CMD)
    CWD = os.getcwd()

    hosts = [Host(host) for host in hosts]

    procs = []
    while len(datasets) > 0:
        ds = datasets.pop(0)
        # load balancing
        hosts.sort()
        host = hosts[0]
        cmd = CMD % ds
        cmd = "ssh %s 'cd %s && %s'" % (host.name, CWD, cmd)
        print "%s: %s" % (host.name, cmd)
        proc = mp.Process(target=run_helper, args=(cmd,))
        proc.start()
        host.njobs += 1
        procs.append(proc)

    for proc in procs:
        proc.join()


if __name__ == "__main__":

    hosts = get_hosts('hosts.sfu.txt')
    hosts = [Host(host) for host in hosts]

    for i in xrange(50):
        hosts.sort()
        hosts[0].njobs += 1
        print ' '.join(map(str, hosts))

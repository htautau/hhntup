import socket
import os
import subprocess
from subprocess import call
import multiprocessing as mp
from higgstautau.datasets import Database
from systematics import iter_systematics


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

    # normalize by the number of CPUs
    cmd = 'python -c "import os; print (os.getloadavg()[0] / open(\\"/proc/cpuinfo\\").read().count(\\"processor\\t:\\"))"'
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


def qsub(cmd, queue='medium', ncpus=1, dry_run=False):

    cmd = "echo '%s' | qsub -q %s -l ncpus=%d" % (
           cmd, queue, ncpus)
    print cmd
    if not dry_run:
        call(cmd, shell=True)


def run(student,
        db,
        datasets,
        hosts,
        nproc=1,
        nice=0,
        setup=None,
        args=None,
        student_args=None,
        use_qsub=False,
        qsub_queue='medium',
        dry_run=False):

    if args is None:
        args = ' '
    else:
        args = ' '.join(args) + ' '

    database = Database(db)

    CMD = "./run -s %s -n %%d --db %s --nice %d %s%%s" % (
            student, db, nice, args)
    if setup is not None:
        CMD = "%s && %s" % (setup, CMD)
    CWD = os.getcwd()

    hosts = [Host(host) for host in hosts]
    datasets = datasets[:]

    procs = []
    while len(datasets) > 0:
        ds = datasets.pop(0)
        # determine actual number of required CPU cores
        files = database[ds].files
        nproc_actual = min(nproc, len(files))
        if not use_qsub:
            # load balancing
            hosts.sort()
            host = hosts[0]
        cmd = CMD % (nproc_actual, ds)
        if student_args is not None:
            cmd = '%s %s' % (cmd, ' '.join(student_args))
        cmd = "cd %s && %s" % (CWD, cmd)
        if use_qsub:
            qsub(cmd, queue=qsub_queue, ncpus=nproc_actual,
                 dry_run=dry_run)
        else: # ssh
            cmd = "ssh %s '%s'" % (host.name, cmd)
            print "%s: %s" % (host.name, cmd)
            if not dry_run:
                proc = mp.Process(target=run_helper, args=(cmd,))
                proc.start()
                procs.append(proc)
            host.njobs += 1

    if not use_qsub and not dry_run:
        for proc in procs:
            proc.join()


def run_systematics(channel, *args, **kwargs):

    for sys_object, sys_type, variation, sys_term in iter_systematics(channel):
        print
        print '======== Running %s systematics ========' % sys_term
        print
        suffix = '--suffix %s' % sys_term
        syst = '--syst-type Systematics.%s --syst-term Systematics.%s.%s' % (
                sys_object, sys_object, sys_term)
        run(*args,
            args=suffix.split(),
            student_args=syst.split(),
            **kwargs)


if __name__ == "__main__":

    hosts = get_hosts('hosts.sfu.txt')
    hosts = [Host(host) for host in hosts]

    while True:
        hosts.sort()
        hosts[0].njobs += 1
        print ' '.join(map(str, hosts))
        print '-' * 10

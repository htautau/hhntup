import time
import socket
import os, errno
import subprocess
from subprocess import call
import multiprocessing as mp
from higgstautau.datasets import Database
from higgstautau import samples
from systematics import iter_systematics
from pbs import qsub
from rootpy.utils.path import mkdir_p


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
    cmd = ('python -c "import os; print (os.getloadavg()[0] / '
           'open(\\"/proc/cpuinfo\\").read().count(\\"processor\\t:\\"))"')
    if not HOSTNAME.startswith(host):
        cmd = "ssh %s '%s'" % (host, cmd)
    load = float(subprocess.Popen(
        cmd, stdout=subprocess.PIPE, shell=True).communicate()[0].strip())
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
        output_path='.',
        setup=None,
        student_args=None,
        use_qsub=False,
        qsub_queue='medium',
        qsub_name_suffix=None,
        dry_run=False,
        separate_student_output=False,
        warnings_as_errors=False,
        **kwargs):

    if not kwargs:
        args = ''
    else:
        args = ' '.join(['--%s %s' % (key, value)
            for key, value in kwargs.items() if value is not None]) + ' '

    if qsub_name_suffix is None:
        qsub_name_suffix = ''
    elif not qsub_name_suffix.startswith('_'):
        qsub_name_suffix = '_' + qsub_name_suffix

    database = Database(db)

    output_path = os.path.normpath(output_path)
    if separate_student_output and os.path.basename(output_path) != student:
        output_path = os.path.join(output_path, os.path.splitext(student)[0])
    if not os.path.exists(output_path):
        if dry_run:
            print "mkdir -p %s" % output_path
        else:
            mkdir_p(output_path)

    python_flags = ''
    if warnings_as_errors:
        python_flags = '-W error'

    CMD = "python %s run --output-path %s -s %s -n %%d --db %s --nice %d %s%%s" % (
           python_flags, output_path, student, db, nice, args)
    if setup is not None:
        CMD = "%s && %s" % (setup, CMD)
    CWD = os.getcwd()

    datasets = datasets[:]

    proc_cmds = []
    while len(datasets) > 0:
        ds = datasets.pop(0)

        output_name = os.path.splitext(student)[0] + '.' + ds
        suffix = kwargs.get('suffix', None)
        if suffix:
            output_name += '_%s' % suffix
        output_name += '.root'
        output_name = os.path.join(output_path, output_name)
        if os.path.exists(output_name):
            print "Output %s already exists. Please delete it and resubmit." % (
                output_name)
            continue

        try:
            files = database[ds].files
        except KeyError:
            print "dataset %s not in database" % ds
            continue

        # determine actual number of required CPU cores
        nproc_actual = min(nproc, len(files))
        cmd = CMD % (nproc_actual, ds)
        if student_args:
            cmd = '%s %s' % (cmd, ' '.join(student_args))
        cmd = "cd %s && %s" % (CWD, cmd)

        if use_qsub: # use the batch system
            qsub(cmd,
                 queue=qsub_queue,
                 ppn=nproc_actual,
                 name=student.strip('.py') + '.' + ds + qsub_name_suffix,
                 stderr_path=output_path,
                 stdout_path=output_path,
                 dry_run=dry_run)

        else: # use simple ssh
            print cmd
            if not dry_run:
                proc_cmds.append(cmd)

    if not use_qsub and not dry_run:
        # use simple ssh with basic load balancing
        hosts = [Host(host) for host in hosts]
        procs = []
        while True:
            active = mp.active_children()
            while len(active) < (2 * len(hosts)) and len(proc_cmds) > 0:
                hosts.sort()
                host = hosts[0]
                cmd = "ssh %s '%s'" % (host.name, proc_cmds.pop(0))
                proc = mp.Process(target=run_helper, args=(cmd,))
                proc.start()
                procs.append(proc)
                host.njobs += 1
                # active_children() joins finished procs
                active = mp.active_children()
            #print time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime())
            #print "jobs: %i running %i queued"%(len(active),len(train_processes))
            if len(proc_cmds) == 0 and len(active) == 0:
                break
            time.sleep(10)


def run_systematics(channel, student, systematics=None, **kwargs):
    student_args = kwargs.pop('student_args', [])
    for sys_variations in iter_systematics(channel):
        if systematics is not None:
            if sys_variations not in systematics:
                continue
        print
        print '======== Running %s systematics ========' % '+'.join(sys_variations)
        print
        syst = '--syst-terms %s' % ','.join(sys_variations)
        run(student,
            student_args=syst.split() + student_args,
            qsub_name_suffix='_'.join(sys_variations),
            suffix='_'.join(sys_variations),
            **kwargs)


def run_systematics_new(channel, student, datasets, systematics,
                        filter_systematics=None, **kwargs):
    student_args = kwargs.pop('student_args', [])
    for sys_variations in systematics:
        if filter_systematics is not None:
            if sys_variations not in filter_systematics:
                continue
        syst = '--syst-terms %s' % ','.join(sys_variations)
        run(student,
            datasets=datasets,
            student_args=syst.split() + student_args,
            qsub_name_suffix='_'.join(sys_variations),
            suffix='_'.join(sys_variations),
            **kwargs)


if __name__ == "__main__":

    hosts = get_hosts('hosts.sfu.txt')
    hosts = [Host(host) for host in hosts]

    while True:
        hosts.sort()
        hosts[0].njobs += 1
        print ' '.join(map(str, hosts))
        print '-' * 10

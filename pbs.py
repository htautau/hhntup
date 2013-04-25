#!/usr/bin/env python
import subprocess
from subprocess import call
import getpass


class Job(object):

    def __init__(self, id, info):

        self.id = id
        self.info = info

    def __getattr__(self, attr):

        return self.info[attr]

    @property
    def name(self):

        return self.info['Job_Name']

    @property
    def health_status(self):

        # is the wall time much higher than the CPU time?
        pass

    @property
    def cputime(self):

        return self.info['resources_used.cput']

    @property
    def walltime(self):

        return self.info['resources_used.walltime']

    @property
    def status(self):

        if 'exec_host' in self.info:
            return self.id, self.info['job_state'], self.info['exec_host'], self.info['Job_Name']
        else:
            return self.id, self.info['job_state'], '-', self.info['Job_Name']


class PBSMonitor(object):

    def __init__(self):

        self.user = getpass.getuser()
        self.jobs = {}
        self.job_names = {}
        self.update()

    def update(self):

        qstat = subprocess.Popen(
            ['qstat', '-f', '-1'],
            stdout=subprocess.PIPE).communicate()[0]
        jobs = qstat.split('\n\n')
        self.jobs = {}
        for block in jobs:
            if not block:
                continue
            block = block.split('\n')
            user = block[2].split(' = ')[-1].split('@')[0]
            if self.user != user:
                continue
            info = {}
            jobid = block[0].split(': ')[-1]
            for line in block[1:]:
                param, value = line.split(' = ')
                info[param.strip()] = value.strip()
            job = Job(jobid, info)
            self.job_names[job.name] = jobid
            self.jobs[jobid] = job

    def has_jobname(self, name):

        return name in self.job_names

    def print_jobs(self):

        for id, job in self.jobs.items():
            print(job.status)


MONITOR = PBSMonitor()


def qsub(cmd,
         queue='medium',
         ncpus=1,
         stderr_path=None,
         stdout_path=None,
         name=None,
         dry_run=False):

    MONITOR.update()
    kwargs = {}
    if name is not None:
        if MONITOR.has_jobname(name):
            print "a job with the name %s already exists" % name
            return
        kwargs['-N'] = name
    if stderr_path is not None:
        kwargs['-e'] = stderr_path
    if stdout_path is not None:
        kwargs['-o'] = stdout_path
    args = ' '.join(['%s %s' % arg for arg in kwargs.items()])
    cmd = "echo '%s' | qsub -q %s %s -l ncpus=%d" % (
           cmd, queue, args, ncpus)
    print cmd
    if not dry_run:
        call(cmd, shell=True)


if __name__ == '__main__':

    MONITOR.print_jobs()

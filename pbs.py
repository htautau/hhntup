import subprocess
from subprocess import call


class PBSMonitor(object):

    def __init__(self):

        self.jobs = {}
        self.update()

    def update(self):

        qstat = subprocess.Popen(
            ['qstat', '-f'],
            stdout=subprocess.PIPE).communicate()[0]
        jobs = qstat.split('\n\n')
        self.jobs = {}
        for block in jobs:
            if not block:
                continue
            block = block.split('\n')
            jobid = block[0][8:]
            jobname = block[1].split('=')[-1].strip()
            self.jobs[jobid] = jobname

    def has_jobname(self, name):

        return name in self.jobs.values()


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


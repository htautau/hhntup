from rootpy.batch import Student, Supervisor

from . import log; log = log[__name__]


class ATLASSupervisor(Supervisor):
    pass


class ATLASStudent(Student):

    def __init__(self, *args, **kwargs):
        super(ATLASStudent, self).__init__(*args, **kwargs)
        self.grl = kwargs.get('grl', None)
        self.events = kwargs.get('events', -1)
        if self.grl:
            log.info("Using GRL: {0}".format(self.grl))

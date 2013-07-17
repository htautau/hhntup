from atlastools.batch import ATLASStudent
from higgstautau import log; log = log[__name__]
from goodruns import GRL


class hhgrl(ATLASStudent):

    def __init__(self, options, **kwargs):

        super(hhgrl, self).__init__(**kwargs)

    def work(self):

        # merge GRL XML strings
        merged_grl = GRL()
        for fname in self.files:
            merged_grl |= GRL('%s:/Lumi/%s' % (fname, self.metadata.treename))
        merged_grl.save('grl.xml')

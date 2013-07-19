from atlastools.batch import ATLASStudent
from rootpy.io import root_open
from higgstautau import log; log = log[__name__]
from goodruns import GRL


class hhgrl(ATLASStudent):

    def __init__(self, options, **kwargs):

        super(hhgrl, self).__init__(**kwargs)

    def work(self):

        # merge GRL XML strings
        grl = GRL()
        for fname in self.files:
            with root_open(fname) as f:
                for key in f.Lumi.keys():
                    grl |= GRL(str(key.ReadObj().GetString()), from_string=True)
        grl.save('grl.xml')

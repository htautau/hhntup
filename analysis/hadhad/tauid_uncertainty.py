#!/usr/bin/env python

from samples import MC_Ztautau
from tauid.p851.selection import selection

"""
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012
Values are percents
"""
EFFIC_UNCERT = {
    'loose': {
        1: 4.,
        3: 8.,
    },
    'medium': {
        1: 5.,
        3: 8.,
    },
    'tight': {
        1: 4.,
        3: 7.,
    },
}


if __name__ == '__main__':

    ztautau = MC_Ztautau()

else:

    def uncertainty(score, pt, prong, nvtx):

        loose = selection('loose', prong, nvtx)
        medium = selection('medium', prong, nvtx)
        tight = selection('tight', prong, nvtx)

        return 0.

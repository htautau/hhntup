#!/usr/bin/env python

from samples import MC_Ztautau
from tauid.p851.selection import selection, nvtx_to_category

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

def uncertainty(score, pt, prong, nvtx):

    loose = selection('loose', prong, nvtx).Eval(pt)
    medium = selection('medium', prong, nvtx).Eval(pt)
    tight = selection('tight', prong, nvtx).Eval(pt)

    if score < loose:
        raise ValueError(
            'No uncertainties defined for scores lower than looose')

    if score < medium:
        return selection_uncertainty('loose', pt, prong, nvtx)
    elif score < tight:
        return selection_uncertainty('medium', pt, prong, nvtx)
    else:
        return selection_uncertainty('tight', pt, prong, nvtx)


def selection_uncertainty(level, pt, prong, nvtx):

    return UNCERT[level][prong][nvtx_to_category(nvtx)].Eval(pt)


if __name__ == '__main__':

    ztautau = MC_Ztautau()

else:

    UNCERT = {}


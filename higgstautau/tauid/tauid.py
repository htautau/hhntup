import os
from rootpy.io import open as ropen
from .selection.p851 import selection, nvtx_to_category, LEVELS, \
    CATEGORIES, PRONGS

HERE = os.path.dirname(os.path.abspath(__file__))

"""
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012
Values are percents / 100
"""
EFFIC_UNCERT_2011 = {
    'loose': {
        1: 0.04,
        3: 0.08,
    },
    'medium': {
        1: 0.05,
        3: 0.08,
    },
    'tight': {
        1: 0.04,
        3: 0.07,
    },
}

"""
https://twiki.cern.ch/twiki/pub/AtlasProtected/TauSystematicsWinterConf2012/tauWG_wtaunu2012winter_01.27.2012_soshi.pdf
"""
EFFIC_SF_2011 = {
    'loose': {
        1: (0.9641, 0.0390),
        3: (0.9415, 0.0755),
    },
    'medium': {
        1: (0.9210, 0.0421),
        3: (0.9849, 0.0800),
    },
    'tight': {
        1: (0.9933, 0.0424),
        3: (0.9711, 0.0710),
    },
}


BDT_UNCERT = {}
with ropen(os.path.join(HERE, 'bdt_uncertainty.root')) as f:
    for level in ('medium', 'tight'):
        BDT_UNCERT[level] = {}
        for prong in PRONGS:
            BDT_UNCERT[level][prong] = {}
            for category in CATEGORIES.keys():
                BDT_UNCERT[level][prong][category] = {}
                for dir in ('high', 'low'):
                    BDT_UNCERT[level][prong][category][dir] = f.Get(
                            '%s_%s_%dp_%s' % (
                                level, dir, prong, category)).Clone()


def uncertainty(score, pt, prong, nvtx):

    if prong > 1:
        prong = 3
    loose = selection('loose', prong, nvtx).Eval(pt)
    medium = selection('medium', prong, nvtx).Eval(pt)
    tight = selection('tight', prong, nvtx).Eval(pt)

    if score < medium:
        raise ValueError(
            'No uncertainties defined for scores lower than medium')

    if score < tight:
        high, low = selection_uncertainty('medium', pt, prong, nvtx)
    else:
        high, low = selection_uncertainty('tight', pt, prong, nvtx)
    return score + high, score - low


def selection_uncertainty(level, pt, prong, nvtx):

    uncert = BDT_UNCERT[level][prong][nvtx_to_category(nvtx)]
    return -1 * uncert['high'].Eval(pt), uncert['low'].Eval(pt)

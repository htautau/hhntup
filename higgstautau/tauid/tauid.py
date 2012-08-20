import os
from rootpy.io import open as ropen
from .p851 import selection, nvtx_to_category, \
    CATEGORIES
from .common import LEVELS, PRONGS, nprong

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


# uncertainty currently only valid for 2011 MC
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

    prong = nprong(prong)
    loose = selection('loose', prong, nvtx).Eval(pt)
    medium = selection('medium', prong, nvtx).Eval(pt)
    tight = selection('tight', prong, nvtx).Eval(pt)

    if score < loose:
        raise ValueError(
            'No uncertainties defined for scores lower than loose')

    high_medium, low_medium = selection_uncertainty('medium', pt, prong, nvtx)
    high_tight, low_tight = selection_uncertainty('tight', pt, prong, nvtx)

    b_m_high = 1. - medium - high_medium
    b_t_high = 1. - tight - high_tight

    b_m_low = 1. - medium + low_medium
    b_t_low = 1. - tight + low_tight

    dx_high = max([high_medium / b_m_high, high_tight / b_t_high]) * (1. - score)
    dx_low = max([low_medium / b_m_low, low_tight / b_t_low]) * (1. - score)

    return score + dx_high, score - dx_low


def selection_uncertainty(level, pt, prong, nvtx):

    uncert = BDT_UNCERT[level][prong][nvtx_to_category(nvtx)]
    return -1 * uncert['high'].Eval(pt), uncert['low'].Eval(pt)

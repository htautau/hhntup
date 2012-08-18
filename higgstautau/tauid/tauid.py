from .selection.p851 import selection, nvtx_to_category, LEVELS, \
    CATEGORIES, PRONGS


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
    }
    'medium': {
        1: (0.9210, 0.0421),
        3: (0.9849, 0.0800),
    }
    'tight': {
        1: (0.9933, 0.0424),
        3: (0.9711, 0.0710),
    }
}


def uncertainty(score, pt, prong, nvtx):

    loose = selection('loose', prong, nvtx).Eval(pt)
    medium = selection('medium', prong, nvtx).Eval(pt)
    tight = selection('tight', prong, nvtx).Eval(pt)

    if score < loose:
        raise ValueError(
            'No uncertainties defined for scores lower than loose')

    if score < medium:
        return selection_uncertainty('loose', pt, prong, nvtx)
    elif score < tight:
        return selection_uncertainty('medium', pt, prong, nvtx)
    else:
        return selection_uncertainty('tight', pt, prong, nvtx)


def selection_uncertainty(level, pt, prong, nvtx):

    return UNCERT[level][prong][nvtx_to_category(nvtx)].Eval(pt)



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

#!/usr/bin/env python

"""
This module defines the separate data periods (analysis-specific)
and associated trigger requirements in each period and
integrated luminosity
"""

from rootpy.extern.tabulartext import PrettyTable


"""
PERIODS = {
    'B-I':  {'runs': (177986, 186493),
             'lumi': 1480.27,
             'triggers': 'EF_tau29_medium1_tau20_medium1'},
    'J':    {'runs': (186516, 186755),
             'lumi': 226.392,
             'triggers': 'EF_tau29_medium1_tau20_medium1'},
    'K':    {'runs': (186873, 187815),
             'lumi': 391.09 + 170.589,
             'triggers': 'EF_tau29_medium1_tau20_medium1'},
    'L-M':  {'runs': (188902, 191933),
             'lumi': 2392.85,
             'triggers': 'EF_tau29_medium1_tau20_medium1 && L1_2TAU11_TAU15'},
}
"""

PERIODS = {
    'B2-M': {'runs': (178044, 191933),
             'lumi': 4604.52}, #4617.}, #4661.2},
}

LUMI = { # pb-1 atlas-lumicalc.cern.ch
    2011: 4604.52, # B2-M skim v9
    2012: 14220.7, #14130.8, # A3-E no run 200841 skim v9
}


if __name__ == '__main__':

    total_lumi = 0
    lumi_table = PrettyTable(['Period(s)', 'Runs', 'Lumi [pb-1]'])
    #trigger_table = PrettyTable(['Period(s)', 'Triggers'])
    for period in sorted(PERIODS.keys()):
        info = PERIODS[period]
        lumi_table.add_row([period, '%i-%i' % info['runs'], '%.3f' % info['lumi']])
        #trigger_table.add_row([period, info['triggers'].replace('||', 'OR').replace('&&', 'AND')])
        total_lumi += info['lumi']
    print lumi_table.get_string(hrules=1)
    print "total lumi: %.3f pb-1" % total_lumi
    #print
    #print trigger_table.get_string(hrules=1)

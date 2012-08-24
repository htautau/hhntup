#!/usr/bin/env python

from higgstautau import tauid
from higgstautau.tauid.p851 import selection, nvtx_to_category, CATEGORIES
from higgstautau.tauid.common import LEVELS, PRONGS
from higgstautau.tauid import EFFIC_UNCERT_2011 as EFFIC_UNCERT

from samples import MC_TauID

import numpy as np
from matplotlib import pyplot as plt
from tauid_uncertainty import efficiency, efficiency_uncertainty


pt = 50000

for prong in tauid.PRONGS:
    loose = selection('loose', prong, 3).Eval(pt)
    medium = selection('medium', prong, 3).Eval(pt)
    tight = selection('tight', prong, 3).Eval(pt)
    scores = np.linspace(0.5, 1., 1000, endpoint=True)
    high = []
    low = []
    for score in scores:
        high_score, low_score = tauid.uncertainty(score, pt, prong, 3)
        high.append(high_score - score)
        low.append(low_score - score)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill_between(scores, low, high, facecolor='yellow', linewidth=0)
    ax.vlines([loose, medium, tight], [-.5, -.5, -.5], [.5, .5, .5], color='k',
            linestyles='--')
    ax.axhline(0., 0., 1.)
    ax.set_xlim(scores[0], scores[-1])
    ax.set_ylim(min(low)*1.2, max(high)*1.2)

    # working point labels
    for wp, loc in (('loose', loose), ('medium', medium), ('tight', tight)):
        plt.text(loc+0.005, ax.get_ylim()[1] / 2., wp,
            ha='left', va='center', color='black', rotation=90)
    plt.ylabel('BDT Score Uncertainty')
    plt.xlabel('BDT Score')
    plt.savefig('tauid_uncertainty_%dp.png' % prong)


# closure test
ztautau = MC_TauID(systematics=False)
for prong in PRONGS:
    for cat_str, category in CATEGORIES.items():
        for level in LEVELS.keys():

            level_selection = selection(level, prong, cat_str)

            uncert_level = EFFIC_UNCERT[level][prong]

            print prong
            print cat_str
            print level
            print "effic: %.6f" % efficiency(ztautau, level_selection, prong, category)
            print "high: %.6f low: %.6f" % efficiency_uncertainty(
                    ztautau, level_selection, prong, category)

            print "=" * 20

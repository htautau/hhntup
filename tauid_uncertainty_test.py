#!/usr/bin/env python

from higgstautau import tauid
from higgstautau.tauid.p851 import selection

import numpy as np
from matplotlib import pyplot as plt


pt = 40000

for prong in tauid.PRONGS:
    loose = selection('loose', prong, 3).Eval(pt)
    medium = selection('medium', prong, 3).Eval(pt)
    tight = selection('tight', prong, 3).Eval(pt)
    scores = np.linspace(loose, 1., 1000, endpoint=True)
    high = []
    low = []
    for score in scores:
        high_score, low_score = tauid.uncertainty(score, pt, prong, 3)
        high.append(high_score - score)
        low.append(low_score - score)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.fill_between(scores, low, high, facecolor='yellow', linewidth=0)
    ax.vlines([medium, tight], [-.5, -.5], [.5, .5], color='k',
            linestyles='dashed')
    ax.set_xlim(loose, 1.)
    ax.set_ylim(min(low)*1.2, max(high)*1.2)
    plt.ylabel('BDT Score Uncertainty')
    plt.xlabel('BDT Score')
    plt.savefig('tauid_uncertainty_%dp.png' % prong)

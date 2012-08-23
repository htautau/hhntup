#!/usr/bin/env python

from higgstautau import tauid
from higgstautau.tauid.p851 import selection
import numpy as np


pt = 40000

for prong in tauid.PRONGS:
    loose = selection('loose', prong, 3).Eval(pt)
    medium = selection('medium', prong, 3).Eval(pt)
    tight = selection('tight', prong, 3).Eval(pt)
    scores = np.linspace(loose, 1., 100, endpoint=True)
    for score in scores:
        print score, tauid.uncertainty(score, pt, prong, 3)

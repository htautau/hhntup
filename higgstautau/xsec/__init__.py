#!/usr/bin/env python

import re
import os

__HERE = os.path.dirname(os.path.abspath(__file__))

SAMPLES = {}

for year, energy in ((11, 7), (12, 8)):
    SAMPLES[year] = {}
    with open(os.path.join(__HERE, '%dTeV.txt' % energy)) as f:
        for line in f.readlines():
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            line = line.split()
            sampleid = int(line[0])
            xsec, kfact, effic = map(float, line[1:-1])
            if sampleid not in SAMPLES[year]:
                SAMPLES[year][sampleid] = {}
            else:
                raise ValueError("duplicate sample {0} in {1}".format(
                    line, f.name))
            SAMPLES[year][sampleid]['lephad'] = {
                'xsec': xsec,
                'effic': effic,
                'kfact': kfact,
                'prod': xsec * kfact / effic,
            }


def xsec_kfact_effic(year, id):
    year %= 1000
    info = SAMPLES[year][id]['lephad']
    return info['xsec'], info['kfact'], info['effic']

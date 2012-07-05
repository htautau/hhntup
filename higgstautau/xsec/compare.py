#!/usr/bin/env python

import re
import os

__HERE = os.path.dirname(os.path.abspath(__file__))

SAMPLE_NAMES = {}

# read sampleid
with open(os.path.join(__HERE, 'sampleid.txt')) as f:
    for line in f.readlines():
        if line.startswith('#'):
            continue
        sampleid, name = line.split()[:2]
        sampleid = int(sampleid)
        SAMPLE_NAMES[sampleid] = name


SAMPLES = {}

# read lephad
with open(os.path.join(__HERE, 'lephad')) as f:
    for line in f.readlines():
        line = line.strip()
        if not line:
            continue
        line = line.split()
        if line[0].startswith('#'):
            if 'kfac' in line[-1]:
                kfact = float(line[-1].split('=')[1].strip(')'))
            else:
                kfact = 1.
        else:
            sampleid, info = line[:2]
            sampleid = int(sampleid)
            info = info.split('*')
            if len(info) >= 2:
                xsec, effic = map(float, info[:2])
            else:
                xsec = float(info[0])
                effic = 1.
            xsec *= 1E3
            if sampleid not in SAMPLES:
                SAMPLES[sampleid] = {}
            SAMPLES[sampleid]['lephad'] = {
                'xsec': xsec,
                'effic': effic,
                'kfact': kfact,
                'prod': xsec * kfact / effic,
            }

# read topmc
with open(os.path.join(__HERE, 'TopMC')) as f:
    for line in f.readlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        sampleid, xsec, kfact = line.split()[:3]
        sampleid = int(sampleid)
        xsec = float(xsec)
        kfact = float(kfact)
        if sampleid not in SAMPLES:
            SAMPLES[sampleid] = {}
        SAMPLES[sampleid]['topmc'] = {
            'xsec': xsec,
            'kfact': kfact,
            'effic': 1.,
            'prod': xsec * kfact,
        }


# read smwz
with open(os.path.join(__HERE, 'SMWZ')) as f:
    for line in f.readlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        sampleid, sample, kfact, effic, xsec = line.split()[:5]
        sampleid = int(sampleid)
        xsec = float(xsec)
        kfact = float(kfact)
        effic =  float(effic)
        if sampleid not in SAMPLES:
            SAMPLES[sampleid] = {}
        SAMPLES[sampleid]['smwz'] = {
            'xsec': xsec,
            'kfact': kfact,
            'effic': effic,
            'prod': xsec * kfact / effic,
        }


def xsec_kfact_effic(group, id):

    info = SAMPLES[id][group]
    return info['xsec'], info['kfact'], info['effic']


if __name__ == '__main__':

    format = '%10s: xs: %-10s k: %-5s eff: %-5s product: %.4g'

    for sampleid in sorted(SAMPLES.keys()):
        if sampleid not in SAMPLE_NAMES:
            continue
        print sampleid, SAMPLE_NAMES[sampleid]
        group_compare = SAMPLES[sampleid]
        for group in sorted(group_compare.keys()):
            group_info = group_compare[group]
            xsec = "%.4g" % group_info['xsec']
            kfact = "%.3f" % group_info['kfact']
            effic = "%.3f" % group_info['effic']
            product = group_info['prod']
            print format % (group, xsec, kfact, effic, product)

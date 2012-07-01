#!/usr/bin/env python

from higgstautau.cutflow import get_parser, make_cutflow


parser = get_parser()
parser.add_argument('--dir', default='../ntuples/midpt')
args = parser.parse_args()

import samples_db

samples = []

for background, info in samples_db.BACKGROUNDS.items():
    for sample in info['samples']:
        sample_info = samples_db.SAMPLES[sample]
        latex = sample_info['latex']
        name = sample_info['name']
        samples.append((latex, name, sample))

make_cutflow(samples, args)

#!/usr/bin/env python

from higgstautau.cutflow import get_parser, make_cutflow


parser = get_parser()
parser.add_argument('--dir', default='../ntuples/midpt')
args = parser.parse_args()

samples = [
    ('B', 'B', 'data-B'),
    ('D', 'D', 'data-D'),
    ('E', 'E', 'data-E'),
    ('F', 'F', 'data-F'),
    ('G', 'G', 'data-G'),
    ('H', 'H', 'data-H'),
    ('I', 'I', 'data-I'),
    ('J', 'J', 'data-J'),
    ('K', 'K', 'data-K'),
    ('L', 'L', 'data-L'),
    ('M', 'M', 'data-M'),
    ('Total', 'Total', 'data-[BDEFGHIJKLM]'),
]

make_cutflow(samples, args)

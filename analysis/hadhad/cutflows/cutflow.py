#!/usr/bin/env python

from higgstautau.cutflow import get_parser, make_cutflow


parser = get_parser()
parser.add_argument('--dir', default='../ntuples/midpt')
parser.add_argument('--short', action='store_true', default=False)
parser.add_argument('--mass', type=int, default=125)
args = parser.parse_args()

samples = [
    ('Data', 'Data', '^data$'),
    (r'$Z\rightarrow\tau_{h}\tau_{h}$', 'Ztautau', 'AlpgenJimmyZtautauNp[\d]_pt20.mc11c'),
    (r'$Z\rightarrow\mu\mu$', 'Zmumu', 'AlpgenJimmyZmumuNp[\d]_pt20.mc11c'),
    (r'$Z\rightarrow ee$', 'Zee', 'AlpgenJimmyZeeNp[\d]_pt20.mc11c'),
    (r'$W\rightarrow\tau_{h}\nu$', 'Wtaunu', 'AlpgenJimmyWtaunuNp[\d]_pt20.mc11c'),
    (r'$W\rightarrow\mu\nu$', 'Wmunu', 'AlpgenJimmyWmunuNp[\d]_pt20.mc11c'),
    (r'$W\rightarrow e\nu$', 'Wenu', 'AlpgenJimmyWenuNp[\d]_pt20.mc11c'),
    (r'ggF $H\rightarrow\tau_{h}\tau_{h}$', 'ggF H tautau', 'PowHegPythia_ggH%d_tautauhh.mc11c' % args.mass),
    (r'VBF $H\rightarrow\tau_{h}\tau_{h}$', 'VBF H tautau', 'PowHegPythia_VBFH%d_tautauhh.mc11c' % args.mass),
]

extra_columns = [2, 3, 5, 6]

if args.short:
    samples = [i for j, i in enumerate(samples) if j not in extra_columns]
    #data_num_format = "%4.3g"
    #num_format = "%4.3g"

make_cutflow(samples, args)

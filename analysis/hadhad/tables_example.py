#!/usr/bin/env python

import tables

f = tables.openFile('../../ntuples/hadhad/HHProcessor/HHProcessor.h5')
t = f.root.data_JetTauEtmiss

#arr = np.array([row[:] for row in f.root.higgstautauhh.where('MET_x > 0')])

arr = t.readWhere('MET_x > 0')

print arr

#ex = tables.Expr('MET_x * 3', uservars=t.colinstances)

# see the readWhere method
#numexpr.evaluate('MET_x * 2', local_dict=dict([(name, arr[name]) for name in
#    arr.dtype.names]))



import tables


f = tables.openFile('asd.h5')

t = tables.root.higgstautauhh


np.array([row[:] for row in f.root.higgstautauhh.where('MET_x > 0')])

ex = tables.Expr('MET_x * 3', uservars=t.colinstances)

#!/usr/bin/env python

import os
from glob import glob
import shutil

DEST = './images/variables'

cats = ['01jet', '2jet']

for cat in cats:
    with open('variables_%s.rst' % cat, 'w') as f:

        plots = glob('../var*%s.png' % cat)
        plots.sort()
        f.write('\n')
        for plot in plots:
            dest = os.path.join(DEST, cat, os.path.basename(plot))
            shutil.copy(plot, dest)
            f.write('.. image:: %s\n' % dest)
            f.write('   :width: 400px\n\n')


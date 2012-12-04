#!/usr/bin/env python

from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-s', '--student', required=True)
parser.add_argument('path')
args = parser.parse_args()

import re
import os
from glob import glob
from rootpy.io import open as ropen
from rootpy.tree import Tree

ds_pattern = re.compile('%s\.(?P<name>.+)\.root$' % (args.student))

with ropen(os.path.join(args.path, args.student) + '.root',
        'RECREATE') as outfile:
    for filename in glob(os.path.join(args.path, '*.root')):
        match = re.search(ds_pattern, filename)
        if not match:
            print "%s is not a valid filename" % filename
            continue
        # replace . and - with _ for natural naming in PyTables
        name = match.group('name').replace('.', '_')
        name = name.replace('-', '_')
        print "Merging in %s ..." % filename
        with ropen(filename, 'READ') as infile:
            intree = infile.higgstautauhh
            outfile.cd()
            outtree = intree.CloneTree(-1, "fast SortBasketsByEntry")
            outtree.OptimizeBaskets()
            outtree.SetName(name)
            outtree.Write()
            # will need to be upated for the next skim
            outcutflow = infile.cutflow.Clone(name=name + '_cutflow')
            outcutflow_event = infile.cutflow_event.Clone(name=name + '_cutflow_event')
            outfile.cd()
            outcutflow.Write()
            outcutflow_event.Write()

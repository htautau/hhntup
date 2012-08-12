#!/usr/bin/env python

from rootpy.io import open as ropen
from rootpy.tree import Cut
from rootpy.plotting import Graph

categories = {
    '_3': Cut('tau_numberOfVertices<=3'),
    '3_5': Cut('3<tau_numberOfVertices<=5'),
    '5_7': Cut('5<tau_numberOfVertices<=7'),
    '7_': Cut('tau_numberOfVertices>7'),
}

levels = {
    1: 'loose',
    2: 'medium',
    3: 'tight'
}

with ropen('bdt_selection.root', 'recreate') as f:

    for prong in (1, 3):
        for cat_str, category in categories.items():
            for level, level_name in levels.items():
                fname = 'sig-bits-%dp-%s-perfB--%d.txt' % (
                        prong, category.safe(parentheses=False), level)
                with open(fname) as fin:
                    lines = fin.readlines()[1:]
                    graph = Graph(len(lines), name='%s_%dp_%s' % (
                        level_name, prong, cat_str))
                    for i, line in enumerate(lines):
                        pt, bdt = map(float, line.strip().split())
                        graph[i] = (pt, bdt)
                    graph.Write()

#!/usr/bin/env python

from rootpy.io import open as ropen
from rootpy.tree import Cut
from rootpy.plotting import Graph
import os

HERE = os.path.dirname(os.path.abspath(__file__))


categories = {
    '_3': Cut('tau_numberOfVertices<=3'),
    '3_5': Cut('3<tau_numberOfVertices<=5'),
    '5_7': Cut('5<tau_numberOfVertices<=7'),
    '7_': Cut('tau_numberOfVertices>7'),
}

levels = {
    'loose': 1,
    'medium': 2,
    'tight': 3,
}


if __name__ == '__main__':

    with ropen(os.path.join(HERE, 'bdt_selection.root'), 'recreate') as f:

        for prong in (1, 3):
            for cat_str, category in categories.items():
                for level_name, level in levels.items():
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
else:

    f = ropen(os.path.join(HERE, 'bdt_selection.root'))


    def nvtx_to_category(nvtx):

        if nvtx <= 3:
            category = '_3'
        elif nvtx <= 5:
            category = '3_5'
        elif nvtx <= 7:
            category = '5_7'
        else:
            category = '7_'
        return category


    def selection(level, prong, nvtx):

        return f.Get('%s_%dp_%s' % (
            level, prong, nvtx_to_category(nvtx)))

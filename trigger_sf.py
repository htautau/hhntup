import os
import ROOT
ROOT.gROOT.SetBatch(True)
from externaltools.bundle_2011 import TauTriggerCorrections
base = TauTriggerCorrections.RESOURCE_PATH

from rootpy.plotting import Canvas, Graph, Legend
from rootpy.plotting.shapes import Line
from rootpy.plotting import get_style, set_style
from rootpy.interactive import wait

import numpy as np

style = get_style('ATLAS', shape='square')
style.SetPadLeftMargin(0.16)
style.SetTitleYOffset(1.6)
style.SetHistTopMargin(0.)
set_style(style)

triggers = [
    ('EF_tau20_medium1', 'red'),
    ('EF_tau20T_medium1', 'black'),
    ('EF_tau29_medium1', 'blue'),
    ('EF_tau29T_medium1', 'green'),
]

pt = np.linspace(20000, 100000, 100)

def draw_curve(_func, name, title, ylow, yhigh):
    graphs = []
    for (trigger, color) in triggers:
        tool = ROOT.TauTriggerCorrections(os.path.join(base, 'triggerSF_%s.root' % trigger))
        func = getattr(tool, _func)
        eff = map(lambda x: func(x, 0), pt)
        eff_low = map(lambda x: func(x, -1), pt)
        eff_high = map(lambda x: func(x, 1), pt)
        graph = Graph(100, name=trigger)
        for i, (p, e, e_low, e_high) in enumerate(zip(pt, eff, eff_low, eff_high)):
            graph.SetPoint(i, p / 1000, e)
            graph.SetPointError(i, 0.4, 0.4, e - e_low, e_high - e)
        graph.linecolor = color
        graph.linewidth = 2
        graph.fillstyle = '/'
        graph.fillcolor = color
        graphs.append(graph)
    c = Canvas()
    leg = Legend(len(graphs),
        pad=c, topmargin=0.4, leftmargin=0.3, textsize=25, margin=0.2)
    for i, g in enumerate(graphs):
        if i == 0:
            g.Draw('3AC')
            g.xaxis.title = 'p_{T} [GeV]'
            g.xaxis.SetLimits(20, 100)
            g.yaxis.SetLimits(ylow, yhigh)
            g.yaxis.SetRangeUser(ylow, yhigh)
            g.yaxis.title = title
        else:
            g.Draw('3C SAME')
        leg.AddEntry(g, g.name, 'L')
    leg.Draw()
    lines = []
    for thresh in (25, 35):
        line = Line(thresh, ylow, thresh, yhigh)
        line.linestyle = 'dashed'
        line.linewidth = 2
        line.Draw()
        lines.append(line)
    c.SaveAs('trigger_{0}.png'.format(name))

draw_curve('getDataEff', 'data_eff', 'Data Efficiency', 0, 1)
draw_curve('getMCEff', 'mc_eff', 'MC Efficiency', 0, 1)
draw_curve('getSF', 'sf', 'Scale Factor', 0.5, 1.3)

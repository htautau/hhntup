import os
import ROOT
ROOT.gROOT.SetBatch(True)
from externaltools import TauTriggerCorrections
base = TauTriggerCorrections.RESOURCE_PATH

from rootpy.plotting import Canvas, Graph, Legend, Hist2D
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

def draw_curve_7(_func, name, title, ylow, yhigh):
    """ Draw 7TeV trigger efficiency curves """
    graphs = []
    for (trigger, color) in triggers:
        tool = ROOT.TauTriggerCorrections(os.path.join(base, 'triggerSF_%s.root' % trigger))
        func = getattr(tool, _func)
        eff = map(lambda x: func(x, 0), pt)
        eff_low = map(lambda x: func(x, -1), pt)
        eff_high = map(lambda x: func(x, 1), pt)
        graph = Graph(len(pt), name=trigger)
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
            g.xaxis.title = '#font[52]{p}_{T} [GeV]'
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
    c.SaveAs('trigger_{0}.eps'.format(name))

# draw 7TeV trigger efficiency curves
draw_curve_7('getDataEff', 'data_eff_7', 'Data Efficiency', 0, 1)
draw_curve_7('getMCEff', 'mc_eff_7', 'MC Efficiency', 0, 1)
draw_curve_7('getSF', 'sf_7', 'MC Efficiency Correction Factor', 0.5, 1.3)

# draw 8TeV trigger efficiency curves
from externaltools import TrigTauEfficiency
base = os.path.join(TrigTauEfficiency.RESOURCE_PATH, 'benchmark_menu')

eveto = 'EVnone'
prong = '1p'
wpflag = 'BDTm'
eta = 0
period = 'periodBD'

triggers = [
    ('EF_tau20Ti_medium1', 'black'),
    ('EF_tau29Ti_medium1', 'red'),
]

pt = np.linspace(20000, 100000, 1000)


def draw_curve_8(_func, name, title, ylow, yhigh, num_errs=1):
    """ Draw 8TeV trigger efficiency curves """
    graphs = []
    for (trigger, color) in triggers:
        tool = ROOT.TrigTauEfficiency()
        tool.loadInputFile(os.path.join(base, 'triggerSF_{0}.root'.format(trigger)))
        func = getattr(tool, _func)
        eff = np.array(map(lambda x: func(x, eta, 0, period, prong, wpflag, eveto), pt))
        errs_low = []
        errs_high = []
        for ierr in xrange(num_errs):
            eff_low = np.array(map(lambda x: func(x, eta, -1, period, prong, wpflag, eveto), pt))
            eff_high = np.array(map(lambda x: func(x, eta, 1, period, prong, wpflag, eveto), pt))
            errs_low.append(eff_low)
            errs_high.append(eff_high)
        # quadrature sum of error
        eff_low = np.sqrt(np.sum([np.power(err, 2) for err in errs_low], axis=0))
        eff_high = np.sqrt(np.sum([np.power(err, 2) for err in errs_high], axis=0))
        graph = Graph(len(pt), name=trigger)
        for i, (p, e, e_low, e_high) in enumerate(zip(pt, eff, eff_low, eff_high)):
            graph.SetPoint(i, p / 1000, e)
            graph.SetPointError(i, 0.4, 0.4, e_low, e_high)
        graph.linecolor = color
        graph.linewidth = 2
        graph.fillstyle = '/'
        graph.fillcolor = color
        graphs.append(graph)
    c = Canvas()
    leg = Legend(len(graphs),
        pad=c, topmargin=0.6, leftmargin=0.3,
        textsize=25, margin=0.2)
    for i, g in enumerate(graphs):
        if i == 0:
            g.Draw('3AL')
            g.xaxis.title = '#font[52]{p}_{T} [GeV]'
            g.xaxis.SetLimits(20, 100)
            g.yaxis.SetLimits(ylow, yhigh)
            g.yaxis.SetRangeUser(ylow, yhigh)
            g.yaxis.title = title
        else:
            g.Draw('3L SAME')
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
    c.SaveAs('trigger_{0}.eps'.format(name))

# draw 8TeV trigger efficiency curves
draw_curve_8('getDataEff', 'data_eff_8', 'Data Efficiency', 0, 1, 2)
draw_curve_8('getMCEff', 'mc_eff_8', 'MC Efficiency', 0, 1, 1)
draw_curve_8('getSF', 'sf_8', 'MC Efficiency Correction Factor', 0.5, 1.3, 3)

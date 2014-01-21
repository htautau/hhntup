import os
import ROOT
from externaltools.bundle_2011 import TauTriggerCorrections
base = TauTriggerCorrections.RESOURCE_PATH

from rootpy.plotting import Canvas, Graph, Legend
from rootpy.plotting import set_style
from rootpy.interactive import wait

import numpy as np

set_style('ATLAS')

triggers = [
    ('triggerSF_EF_tau20_medium1', 'red'),
    ('triggerSF_EF_tau29_medium1', 'blue'),
    ('triggerSF_EF_tau20T_medium1', 'black'),
    ('triggerSF_EF_tau29T_medium1', 'green'),
]

pt = np.linspace(20000, 100000, 100)
graphs = []

for (trigger, color) in triggers:
    tool = ROOT.TauTriggerCorrections(os.path.join(base, '%s.root' % trigger))
    eff = map(lambda x: tool.getDataEff(x, 0), pt)
    graph = Graph(100, name=trigger)
    for i, (p, e) in enumerate(zip(pt, eff)):
        graph.SetPoint(i, p / 1000, e)
    graph.markercolor = color
    graphs.append(graph)

c = Canvas()
leg = Legend(len(graphs), pad=c, topmargin=0.4, leftmargin=0.2, textsize=25)
for i, g in enumerate(graphs):
    if i == 0:
        g.Draw('AP')
        g.xaxis.title = 'p_{T} [GeV]'
        g.yaxis.title = 'Data Efficiency'
    else:
        g.Draw('P SAME')
    leg.AddEntry(g, g.name)
leg.Draw()
c.SaveAs('trigger_sf.png')
wait()

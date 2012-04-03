from rootpy.io import open
from rootpy.plotting import Hist, Canvas
from rootpy.tree import Tree
from atlastools import style as  atlastyle
from rootpy.common import logon, set_style
logon(batch=True)
set_style(atlastyle.get_style())

import ROOT
ROOT.gStyle.SetLineWidth(2)

f = open("../../gridcode/VBFH130hh.root")
t = f.Get("VBFH130hh")
t_mc = f.Get("VBFH130hh_mc")


h1 = Hist(100, -6, 6)
h1.SetXTitle("True tau #eta")

h2 = h1.Clone()
h2.SetXTitle("Reco tau #eta")

h3 = Hist(100, -6, 6)
h3.SetXTitle("True tau boosted #eta")

h4 = h1.Clone()
h4.SetXTitle("Reco tau boosted #eta")

t_mc.Draw("fourvect_vis.Eta()","fourvect_vis.Pt()>20000",hist=h1)
t_mc.Draw("fourvect_vis_boosted.Eta()","fourvect_vis.Pt()>20000",hist=h3)

h1.SetFillColor('red')
h1.SetFillStyle('\\')
h3.SetFillColor('red')
h3.SetFillStyle('\\')

c = Canvas()
h1.Draw("hist")

c.SaveAs("truetau_eta.png")

h3.Draw("hist")

c.SaveAs("truetau_eta_boosted.png")

c.Clear()
t.Draw("tau1_fourvect.Eta()","tau1_matched",hist=h2)
t.Draw("tau2_fourvect.Eta()","tau2_matched",hist=h2)

t.Draw("tau1_fourvect_boosted.Eta()","tau1_matched",hist=h4)
t.Draw("tau2_fourvect_boosted.Eta()","tau2_matched",hist=h4)

h2.SetFillColor('red')
h2.SetFillStyle('\\')
h4.SetFillColor('red')
h4.SetFillStyle('\\')

h2.Draw("hist")

c.SaveAs("tau_eta.png")

h4.Draw("hist")

c.SaveAs("tau_eta_boosted.png")

"""
h3 = Hist(100, 0, 1)
h3.SetXTitle("Jet Vertex Fraction")
t.Draw("jet_AntiKt4TopoEM_jvtxf","selected && jet_AntiKt4TopoEM_jvtxf > 0",hist=h3)
c.Clear()
h3.Draw("hist")
c.SaveAs("jvf.eps")
"""

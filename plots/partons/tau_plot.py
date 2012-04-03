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

h3 = h1.Clone()
h3.SetXTitle("Reco tau #eta")

t_mc.Draw("eta_vis","pt_vis>20000",hist=h1)

h1.SetFillColor('red')
h1.SetFillStyle('\\')

c = Canvas()
h1.Draw("hist")

c.SaveAs("truetau_eta.png")

c.Clear()
t.Draw("tau1_seedCalo_eta","tau1_matched && selected",hist=h3)
t.Draw("tau2_seedCalo_eta","tau2_matched && selected",hist=h3)

h3.SetFillColor('red')
h3.SetFillStyle('\\')

h3.Draw("hist")

c.SaveAs("tau_eta.png")


"""
h3 = Hist(100, 0, 1)
h3.SetXTitle("Jet Vertex Fraction")
t.Draw("jet_AntiKt4TopoEM_jvtxf","selected && jet_AntiKt4TopoEM_jvtxf > 0",hist=h3)
c.Clear()
h3.Draw("hist")
c.SaveAs("jvf.eps")
"""

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

f_sm = open("../../gridcode/VBFH130hh_SM.root")
t_sm = f_sm.Get("VBFH130hh_SM")


h1 = Hist(100, -6, 6)
h1.SetXTitle("True quark #eta")
h2 = h1.Clone()
h2.SetXTitle("True quark #eta")

h3 = h1.Clone()
h3.SetXTitle("Reco jet #eta")
h4 = h1.Clone()
h4.SetXTitle("Reco jet #eta")

t_sm.Draw("parton1_eta","parton1_pt > 20000 && selected",hist=h1)
t_sm.Draw("parton2_eta","parton2_pt > 20000 && selected",hist=h2)

t_sm.Draw("jet1_eta","jet1_pt > 20000 && selected",hist=h3)
t_sm.Draw("jet2_eta","jet2_pt > 20000 && selected",hist=h4)

h1.SetFillColor('red')
h1.SetFillStyle('\\')
h2.SetFillColor('blue')
h2.SetFillStyle('/')

c = Canvas()
h1.Draw("hist")
h2.Draw("hist same")

c.SaveAs("parton_eta.png")

c.Clear()
t_sm.Draw("jet1_eta","jet2_eta > -6 && selected",hist=h3)
t_sm.Draw("jet2_eta","jet2_eta > -6 && selected",hist=h4)

h3.SetFillColor('red')
h3.SetFillStyle('\\')
h4.SetFillColor('blue')
h4.SetFillStyle('/')

h3.Draw("hist")
h4.Draw("hist same")

c.SaveAs("jet_eta.png")


"""
h3 = Hist(100, 0, 1)
h3.SetXTitle("Jet Vertex Fraction")
t.Draw("jet_AntiKt4TopoEM_jvtxf","selected && jet_AntiKt4TopoEM_jvtxf > 0",hist=h3)
c.Clear()
h3.Draw("hist")
c.SaveAs("jvf.eps")
"""

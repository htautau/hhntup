#!/usr/bin/env python
from rootpy.io import open
from rootpy.plotting import Hist, Canvas
from rootpy.tree import Tree
from atlastools import style as  atlastyle
from rootpy.common import logon, set_style
logon(batch=True)
set_style(atlastyle.get_style())

import ROOT
ROOT.gStyle.SetLineWidth(2)

f = open("/global/endw/higgs/data/data.root")
t_data = f.Get("data")

f_sm = open("../../../gridcode/VBFH130hh_SM.root")
t_sm = f_sm.Get("VBFH130hh_SM")


h1 = Hist(100, -6, 6)
h1.SetXTitle("Data reco jet #eta")
h2 = h1.Clone()
h2.SetXTitle("Data reco jet #eta")

h3 = h1.Clone()
h3.SetXTitle("Reco jet #eta")
h4 = h1.Clone()
h4.SetXTitle("Reco jet #eta")

h5 = h1.Clone()
h5.SetXTitle("Data reco jet boosted #eta")
h6 = h1.Clone()
h6.SetXTitle("Data reco jet boosted #eta")

h7 = h1.Clone()
h7.SetXTitle("Reco jet boosted #eta")
h8 = h1.Clone()
h8.SetXTitle("Reco jet boosted #eta")

t_data.Draw("jet1_fourvect.Eta()","jet1_fourvect.Pt() > 20000 && tau1_BDTJetScore < 0.5",hist=h1)
t_data.Draw("jet2_fourvect.Eta()","jet2_fourvect.Pt() > 20000 && tau1_BDTJetScore < 0.5",hist=h2)

t_sm.Draw("jet1_fourvect.Eta()","jet1_fourvect.Pt() > 20000",hist=h3)
t_sm.Draw("jet2_fourvect.Eta()","jet2_fourvect.Pt() > 20000",hist=h4)

t_data.Draw("jet1_fourvect_boosted.Eta()","jet1_fourvect.Pt() > 20000 && tau1_BDTJetScore < 0.5",hist=h5)
t_data.Draw("jet2_fourvect_boosted.Eta()","jet2_fourvect.Pt() > 20000 && tau1_BDTJetScore < 0.5",hist=h6)

t_sm.Draw("jet1_fourvect_boosted.Eta()","jet1_fourvect.Pt() > 20000",hist=h7)
t_sm.Draw("jet2_fourvect_boosted.Eta()","jet2_fourvect.Pt() > 20000",hist=h8)

h1.SetFillColor('red')
h1.SetFillStyle('\\')
h2.SetFillColor('blue')
h2.SetFillStyle('/')

c = Canvas()
h1.Draw("hist")
h2.Draw("hist same")

c.SaveAs("parton_eta.png")

h5.SetFillColor('red')
h5.SetFillStyle('\\')
h6.SetFillColor('blue')
h6.SetFillStyle('/')

c = Canvas()
h5.Draw("hist")
h6.Draw("hist same")

c.SaveAs("parton_eta_boosted.png")

c.Clear()

h3.SetFillColor('red')
h3.SetFillStyle('\\')
h4.SetFillColor('blue')
h4.SetFillStyle('/')

h3.Draw("hist")
h4.Draw("hist same")

c.SaveAs("jet_eta.png")

c.Clear()
h7.SetFillColor('red')
h7.SetFillStyle('\\')
h8.SetFillColor('blue')
h8.SetFillStyle('/')

h7.Draw("hist")
h8.Draw("hist same")

c.SaveAs("jet_eta_boosted.png")

"""
h3 = Hist(100, 0, 1)
h3.SetXTitle("Jet Vertex Fraction")
t.Draw("jet_AntiKt4TopoEM_jvtxf","selected && jet_AntiKt4TopoEM_jvtxf > 0",hist=h3)
c.Clear()
h3.Draw("hist")
c.SaveAs("jvf.eps")
"""

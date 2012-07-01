#!/usr/bin/env python

from rootpy.io import open
from rootpy.plotting import Hist, Canvas
from rootpy.tree import Tree
from atlastools import style as  atlastyle
from rootpy.common import logon, set_style, draw
logon(batch=True)
set_style(atlastyle.get_style())

import ROOT
ROOT.gStyle.SetLineWidth(2)

f = open("../../VBFTruthProcessor.PowHegPythia_VBFH130_tautauhh.root")
t = f.Get("PowHegPythia_VBFH130_tautauhh")

"""
Flagging stats
"""
total = t.GetEntries()
print "Total events: %i" % total
print "2 taus + 2 jets flagged: %.3f%%" % (100 * t.GetEntries('jet2_fourvect.Pt()>0') / float(total))

print "2 visible true taus: %.3f%%" % (100 * t.GetEntries('trueTau1_fourvect_vis.Et() > 10000 && abs(trueTau1_fourvect_vis.Pt()) < 2.5 && '
                                                          'trueTau2_fourvect_vis.Et() > 10000 && abs(trueTau2_fourvect_vis.Pt()) < 2.5') / float(total))
print "2 taus truth-matched: %.3f%%" % (100 * t.GetEntries('tau1_matched && tau2_matched') / float(total))
print "2 jets truth-matched: %.3f%%" % (100 * t.GetEntries('jet1_matched && jet2_matched') / float(total))
print "All truth-matched: %.3f%%" % (100 * t.GetEntries('jet1_matched && jet2_matched && tau1_matched && tau2_matched') / float(total))

"""
Jets and true jets
"""

h1 = Hist(100, 10, 300, format='hist', legendstyle='F', title='True quark')

h3 = h1.Clone(title='Reco jet')

h5 = h1.Clone()

h7 = h3.Clone()

t.Draw("parton1_fourvect.Pt()/1000","parton1_fourvect.Pt() > 20000",hist=h1)
t.Draw("parton2_fourvect.Pt()/1000","parton2_fourvect.Pt() > 20000",hist=h1)

t.Draw("jet1_fourvect.Pt()/1000","jet1_fourvect.Pt() > 20000",'error==0',hist=h3)
t.Draw("jet2_fourvect.Pt()/1000","jet2_fourvect.Pt() > 20000",'error==0',hist=h3)

t.Draw("parton1_fourvect_boosted.Pt()/1000","parton1_fourvect.Pt() > 20000",hist=h5)
t.Draw("parton2_fourvect_boosted.Pt()/1000","parton2_fourvect.Pt() > 20000",hist=h5)

t.Draw("jet1_fourvect_boosted.Pt()/1000","jet1_fourvect.Pt() > 20000",'error==0',hist=h7)
t.Draw("jet2_fourvect_boosted.Pt()/1000","jet2_fourvect.Pt() > 20000",'error==0',hist=h7)

h1.SetFillColor('red')
h1.SetFillStyle('\\')

h3.SetFillColor('red')
h3.SetFillStyle('-')

c = Canvas()

draw([h1, h3], axislabels=['p_{T} [GeV]', 'Events'], greedylegend=True)

c.SaveAs("jet_parton_pt.png")

c.Clear()

h5.SetFillColor('red')
h5.SetFillStyle('\\')

h7.SetFillColor('red')
h7.SetFillStyle('-')

c = Canvas()

draw([h5, h7], axislabels=['p_{T} [GeV] boosted', 'Events'], greedylegend=True)

c.SaveAs("jet_parton_pt_boosted.png")


"""
Taus and true taus
"""

h1 = Hist(100, 10, 300, format='hist', legendstyle='F', title='True tau')
h2 = h1.Clone(title='Reco tau')
h3 = h1.Clone()
h4 = h2.Clone()
h5 = h2.Clone()

t.Draw("tau1_fourvect.Pt()/1000",'error==0',hist=h2)
t.Draw("tau2_fourvect.Pt()/1000",'error==0',hist=h2)

t.Draw("trueTau1_fourvect_vis.Pt()/1000","error==0 && trueTau1_fourvect_vis.Pt()>20000",hist=h1)
t.Draw("trueTau2_fourvect_vis.Pt()/1000","error==0 && trueTau2_fourvect_vis.Pt()>20000",hist=h1)

t.Draw("tau1_fourvect_boosted.Pt()/1000",'error==0',hist=h4)
t.Draw("tau2_fourvect_boosted.Pt()/1000",'error==0',hist=h4)

t.Draw("trueTau1_fourvect_vis_boosted.Pt()/1000","error==0 && trueTau1_fourvect_vis.Pt()>20000",hist=h3)
t.Draw("trueTau2_fourvect_vis_boosted.Pt()/1000","error==0 && trueTau2_fourvect_vis.Pt()>20000",hist=h3)

h1.SetFillColor('red')
h1.SetFillStyle('\\')
h2.SetFillColor('blue')
h2.SetFillStyle('/')

c = Canvas()

draw([h1, h2], axislabels=['p_{T} [GeV]', 'Events'], greedylegend=True)

c.SaveAs("tau_pt.png")

h5.SetFillColor('blue')
h5.SetFillStyle('/')

c.Clear()

h3.SetFillColor('red')
h3.SetFillStyle('\\')
h4.SetFillColor('blue')
h4.SetFillStyle('/')

draw([h3, h4], axislabels=['p_{T} [GeV] boosted', 'Events'], greedylegend=True)

c.SaveAs("tau_pt_boosted.png")

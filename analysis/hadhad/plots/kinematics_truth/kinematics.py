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

print "2 visible true taus: %.3f%%" % (100 * t.GetEntries('trueTau1_fourvect_vis.Et() > 10000 && abs(trueTau1_fourvect_vis.Eta()) < 2.5 && '
                                                          'trueTau2_fourvect_vis.Et() > 10000 && abs(trueTau2_fourvect_vis.Eta()) < 2.5') / float(total))
print "2 taus truth-matched: %.3f%%" % (100 * t.GetEntries('tau1_matched && tau2_matched') / float(total))
print "2 jets truth-matched: %.3f%%" % (100 * t.GetEntries('jet1_matched && jet2_matched') / float(total))
print "All truth-matched: %.3f%%" % (100 * t.GetEntries('jet1_matched && jet2_matched && tau1_matched && tau2_matched') / float(total))

"""
Jets and true jets
"""

h1 = Hist(100, -6, 6, format='hist', legendstyle='F', title='True quark')
h2 = h1.Clone()

h3 = h1.Clone(title='Reco jet')
h4 = h1.Clone(title='Reco jet')

h5 = h1.Clone()
h6 = h1.Clone()

h7 = h3.Clone()
h8 = h3.Clone()

t.Draw("parton1_fourvect.Eta()","parton1_fourvect.Pt() > 20000",hist=h1)
t.Draw("parton2_fourvect.Eta()","parton2_fourvect.Pt() > 20000",hist=h2)

t.Draw("jet1_fourvect.Eta()","jet1_fourvect.Pt() > 20000",'error==0',hist=h3)
t.Draw("jet2_fourvect.Eta()","jet2_fourvect.Pt() > 20000",'error==0',hist=h4)

t.Draw("parton1_fourvect_boosted.Eta()","parton1_fourvect.Pt() > 20000",hist=h5)
t.Draw("parton2_fourvect_boosted.Eta()","parton2_fourvect.Pt() > 20000",hist=h6)

t.Draw("jet1_fourvect_boosted.Eta()","jet1_fourvect.Pt() > 20000",'error==0',hist=h7)
t.Draw("jet2_fourvect_boosted.Eta()","jet2_fourvect.Pt() > 20000",'error==0',hist=h8)

h1.SetFillColor('red')
h1.SetFillStyle('\\')
h2.SetFillColor('blue')
h2.SetFillStyle('/')

h3.SetFillColor('red')
h3.SetFillStyle('-')
h4.SetFillColor('blue')
h4.SetFillStyle('|')

c = Canvas()

draw([h1, h2, h3, h4], axislabels=['#eta', 'Events'])

c.SaveAs("jet_parton_eta.png")

c.Clear()

h5.SetFillColor('red')
h5.SetFillStyle('\\')
h6.SetFillColor('blue')
h6.SetFillStyle('/')

h7.SetFillColor('red')
h7.SetFillStyle('-')
h8.SetFillColor('blue')
h8.SetFillStyle('|')

c = Canvas()

draw([h5, h6, h7, h8], axislabels=['#eta boosted', 'Events'])

c.SaveAs("jet_parton_eta_boosted.png")


"""
Taus and true taus
"""

h1 = Hist(100, -6, 6, format='hist', legendstyle='F', title='True tau')
h2 = h1.Clone(title='Reco tau')
h3 = h1.Clone()
h4 = h2.Clone()
h5 = h2.Clone()

t.Draw("tau1_fourvect.Eta()",'error==0',hist=h2)
t.Draw("tau2_fourvect.Eta()",'error==0',hist=h2)

t.Draw("trueTau1_fourvect_vis.Eta()","error==0 && trueTau1_fourvect_vis.Pt()>20000",hist=h1)
t.Draw("trueTau2_fourvect_vis.Eta()","error==0 && trueTau2_fourvect_vis.Pt()>20000",hist=h1)

t.Draw("tau1_fourvect_boosted.Eta()",'error==0',hist=h4)
t.Draw("tau2_fourvect_boosted.Eta()",'error==0',hist=h4)

t.Draw("trueTau1_fourvect_vis_boosted.Eta()","error==0 && trueTau1_fourvect_vis.Pt()>20000",hist=h3)
t.Draw("trueTau2_fourvect_vis_boosted.Eta()","error==0 && trueTau2_fourvect_vis.Pt()>20000",hist=h3)

h1.SetFillColor('red')
h1.SetFillStyle('\\')
h2.SetFillColor('blue')
h2.SetFillStyle('/')

c = Canvas()

draw([h1, h2], axislabels=['#eta', 'Events'])

c.SaveAs("tau_eta.png")

h5.SetFillColor('blue')
h5.SetFillStyle('/')

c.Clear()

h3.SetFillColor('red')
h3.SetFillStyle('\\')
h4.SetFillColor('blue')
h4.SetFillStyle('/')

draw([h3, h4], axislabels=['#eta boosted', 'Events'])

c.SaveAs("tau_eta_boosted.png")


jet_mass = Hist(100, 10, 600, title='true diquark mass', format='hist', fillstyle='/', legendstyle='F')
tau_mass = jet_mass.Clone(title='true ditau visible mass', fillstyle='\\', legendstyle='F')

jet_mass.SetFillStyle('/')
jet_mass.SetFillColor('red')
tau_mass.SetFillStyle('\\')
tau_mass.SetFillColor('blue')

t.Draw("true_Mvis_tau1_tau2/1000", '(trueTau1_fourvect_vis.Pt()>20000) && (trueTau2_fourvect_vis.Pt()>20000) && error == 0',hist=tau_mass)
t.Draw("true_M_quark1_quark2/1000", '(parton1_fourvect.Pt() > 20000) && (parton2_fourvect.Pt() > 20000) && error == 0',hist=jet_mass)

c.Clear()

draw([jet_mass, tau_mass], axislabels=['Visible Mass [GeV]', 'Events'], greedylegend=True)


c.SaveAs('quark_tau_mass.png')

jet_mass = Hist(100, 10, 600, title='reco VBF dijet mass', format='hist', fillstyle='/', legendstyle='F')
tau_mass = jet_mass.Clone(title='reco ditau visible mass', fillstyle='\\', legendstyle='F')

jet_mass.SetFillStyle('/')
jet_mass.SetFillColor('red')
tau_mass.SetFillStyle('\\')
tau_mass.SetFillColor('blue')

t.Draw("Mvis_tau1_tau2/1000", 'error == 0',hist=tau_mass)
t.Draw("M_jet1_jet2/1000", 'error == 0',hist=jet_mass)

c.Clear()

draw([jet_mass, tau_mass], axislabels=['Visible Mass [GeV]', 'Events'], greedylegend=True)


c.SaveAs('reco_jet_tau_mass.png')


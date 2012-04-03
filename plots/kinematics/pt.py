#!/usr/bin/env python

import ROOT
from rootpy.common import logon, makeLabel
from rootpy.plotting import Canvas, Hist
from rootpy.io import openFile
from rootpy.tree import Cut
logon(batch=True)

f = openFile("/global/endw/higgs/data/data.root")
t = f.Get("data")
c = Canvas()
c.cd()
mass = Hist(100,20,200)

from numpy import linspace

cuts = linspace(20000, 100000, 200)

BDTcut = 0.7

for i, pTcut in enumerate(cuts):
    
    print "frame %i"% i
    cut = Cut("(tau1_charge * tau2_charge == -1 ) && tau1_BDTJetScore>%f && tau2_BDTJetScore>%f" % (BDTcut, BDTcut))
    cut &= "tau1_pt>%f && tau2_pt>%f" % (pTcut, pTcut)

    c.Clear()
    mass.Reset()
    t.Draw("Mvis_tau1_tau2/1000", cut, hist=mass)
    mass.SetXTitle("M_{vis}(#tau_{1},#tau_{2}) [GeV]")
    mass.SetYTitle("Events")
    mass.Draw()
    label = makeLabel(0.5,0.8,"BDT > %.2f, p_{T} > %.0f GeV" % (BDTcut, pTcut/1000), size=25)
    label.Draw()
    c.SaveAs("masspt/mass_%04d.png"% i)

    """
    c.Clear()
    tracks.Reset()
    t.Draw("tau1_numTrack", cut, hist=tracks)
    tracks.SetXTitle("Number of Tracks")
    tracks.SetYTitle("#tau Candidates")
    tracks.Draw()
    label = routines.makeLabel(0.7,0.8,"BDT > %.2f"% BDTcut, size=25)
    label.Draw()
    c.SaveAs("tracks_%04d.png"% i)
    """

"""
h = Hist(50,80,300)
t.Draw("Mvis_tau1_tau2/1000", cut, hist=h)
h.SetXTitle("M_{vis}(#tau_{1},#tau_{2}) [GeV]")
h.Draw()
c.SaveAs("massroi.png")

from rootpy import root2matplotlib as r2m
from matplotlib import pyplot as plt

plt.figure()
r2m.errorbar(h)
plt.savefig("massplot.png")
"""

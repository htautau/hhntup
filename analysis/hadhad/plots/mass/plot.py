#!/usr/bin/env python

import ROOT
from rootpy import common
from atlastools import style
from rootpy.plotting import Canvas, Hist
from rootpy.io import open as openFile
from rootpy.tree import Cut
common.logon(batch=True)
common.set_style(style.get_style())


f = openFile("/global/endw/higgs/data/data.root")
t = f.Get("data")
f2 = openFile("/global/endw/higgs/Ztautau/Ztautau.root")
t2 = f2.Get("Ztautau")

c = Canvas()
c.cd()
mass = Hist(100,20,200)
mass2 = Hist(100,20,200)
mass2.SetFillStyle("/")
mass2.SetFillColor("red")

MMC_mass = Hist(100,20,400)

from numpy import logspace
cuts = list(reversed(1 - (logspace(0,1.05,150) - 1)/10))[2:]
cuts.insert(0,0.)

var, labelname = "Mvis_tau1_tau2/1000", 'vis'
var, labelname = "MMC_mass", 'MMC'

for i, BDTcut in enumerate(cuts):
    
    print "frame %i"% i
    cut = Cut("(tau1_charge * tau2_charge == -1 ) && tau1_BDTJetScore>%f && tau2_BDTJetScore>%f" % (BDTcut, BDTcut))

    c.Clear()
    mass.Reset()
    mass2.Reset()
    t.Draw(var, cut, hist=mass)
    t2.Draw(var, cut & "tau1_matched && tau2_matched && EF_tau29_medium1_tau20_medium1", hist=mass2)
    mass.SetXTitle("M_{%s}(#tau_{1},#tau_{2}) [GeV]" % labelname)
    mass.SetYTitle("Events")
    mass /= sum(mass)
    mass2 /= sum(mass2)
    
    _max = 1.2 * max(max(mass), max(mass2))
    mass.SetMinimum(0)
    mass2.SetMinimum(0)
    mass.SetMaximum(_max)
    mass2.SetMaximum(_max)

    mass.Draw()
    mass2.Draw("hist same")
    label = common.makeLabel(0.7,0.8,"BDT > %.2f"% BDTcut, size=25)
    label.Draw()
    c.SaveAs("mass_overlay/mass_%s_%04d.png"% (labelname, i))
    
    """
    c.Clear()
    MMC_mass.Reset()
    t.Draw("MMC_mass", cut, hist=MMC_mass)
    MMC_mass.SetXTitle("M(#tau_{1},#tau_{2}) [GeV]")
    MMC_mass.SetYTitle("Events")
    MMC_mass.Draw()
    label = routines.makeLabel(0.7,0.8,"BDT > %.2f"% BDTcut, size=25)
    label.Draw()
    c.SaveAs("mass/MMC_mass_%04d.eps"% i)
    """

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

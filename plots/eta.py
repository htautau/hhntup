#!/usr/bin/env python

from rootpy.io import openFile
from rootpy.plotting import Canvas, Hist
from rootpy.tree import Tree, Cut
from rootpy.common import logon

logon(batch=True)

mc_file = openFile("/global/endw/higgs/VBFH130hh/VBFH130hh.root")
mc = mc_file.Get("VBFH130hh")
data_file = openFile("/global/endw/higgs/data/data.root")
data = data_file.Get("data")

c = Canvas()

mmc_mc = Hist(100, 0, 10)
mmc_data = Hist(100, 0, 10)
mmc_data.SetXTitle("#Delta#eta(jet_{1}, jet_{2})")
mmc_data.SetYTitle("A.U.")

selection = Cut("tau1_BDTJetScore > .7 && tau2_BDTJetScore > .7")
selection &= "jet1_E>40000 && jet2_E>40000"
mc_selection = Cut("selected && tau1_matched && tau2_matched")

mc.Draw("jetDeltaEta", mc_selection & selection, hist=mmc_mc)
data.Draw("jetDeltaEta", selection, hist=mmc_data)

mmc_mc /= sum(mmc_mc)
mmc_data /= sum(mmc_data)

_max = 1.1*max((max(mmc_mc), max(mmc_data)))
mmc_mc.SetMaximum(_max)
mmc_data.SetMaximum(_max)

mmc_data.Draw()
mmc_mc.Draw("hist same")

c.SaveAs("deta.png")



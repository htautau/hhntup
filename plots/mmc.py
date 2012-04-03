#!/usr/bin/env python

from rootpy.io import open
from rootpy.plotting import Canvas, Hist
from rootpy.tree import Tree, Cut
from rootpy.common import logon

logon(batch=True)

mc_file = open("/global/endw/higgs/VBFH130hh/VBFH130hh.root")
mc = mc_file.Get("VBFH130hh")
data_file = open("/global/endw/higgs/data/data.root")
data = data_file.Get("data")

c = Canvas()

mmc_mc = Hist(300, 0, 300)
mmc_data = Hist(300, 0, 300)
mmc_data_low = Hist(50, 0, 400)
mmc_data.SetXTitle("M_{MMC} [GeV]")
mmc_data.SetYTitle("A.U.")

mmc_selection = Cut("MMC_mass > 0")
selection = Cut("tau1_BDTJetScore > .7 && tau2_BDTJetScore > .7")
mc_selection = Cut("selected && tau1_matched && tau2_matched")

mc.Draw("MMC_mass", mc_selection & selection & mmc_selection, hist=mmc_mc)
data.Draw("MMC_mass", selection & mmc_selection, hist=mmc_data)
data.Draw("MMC_mass", selection & mmc_selection, hist=mmc_data_low)


print sum(mmc_data_low)

mmc_mc /= sum(mmc_mc)
mmc_data /= sum(mmc_data)

_max = 1.1*max((max(mmc_mc), max(mmc_data)))
mmc_mc.SetMaximum(_max)
mmc_data.SetMaximum(_max)

mmc_data.Draw("hist")
mmc_mc.Draw("hist same")

c.SaveAs("mass_mmc.eps")

c.Clear()

m_mc = Hist(300, 0, 200)
m_data = Hist(300, 0, 200)
m_data.SetXTitle("M_{vis} [GeV]")
m_data.SetYTitle("A.U.")

mc.Draw("Mvis_tau1_tau2/1000", mc_selection & selection, hist=m_mc)
data.Draw("Mvis_tau1_tau2/1000", selection, hist=m_data)

m_mc /= sum(m_mc)
m_data /= sum(m_data)

_max = 1.1*max((max(m_mc), max(m_data)))
m_mc.SetMaximum(_max)
m_data.SetMaximum(_max)

m_data.Draw("hist")
m_mc.Draw("hist same")

c.SaveAs("mass_vis.png")

c.Clear()
mmc_data_low.Draw("hist")
c.SaveAs("mass_mmc_low.png")

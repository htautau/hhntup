#!/usr/bin/env python

import ROOT
from rootpy.utils import asrootpy
from rootpy.io import open as ropen
from rootpy.tree import Cut, Tree
from rootpy.plotting import Hist, Canvas
from rootpy import common
from atlastools import style as atlasstyle
import os


common.logon(batch=True)
common.set_style(style = atlasstyle.get_style())

data_file = ropen('../data_M.root')
signal_file = ropen('../PowHegPythia_VBFH130_tautauhh.root')

data = data_file.data_M
signal = signal_file.PowHegPythia_VBFH130_tautauhh

data.SetWeight(1)
signal.SetWeight(1)

canvas = Canvas(800, 600)

output = ropen('VBF_BDT.root', 'RECREATE')
factory = ROOT.TMVA.Factory('VBF_BDT', output, "V:!Silent:!Color:DrawProgressBar")

tmp = ropen('tmp.root', 'RECREATE')

SS = Cut('tau1_charge * tau2_charge == 1')
OS = Cut('tau1_charge * tau2_charge == -1')
loose_taus = Cut('tau1_JetBDTSigLoose && tau2_JetBDTSigLoose')

# opposite sign
signal_selected = asrootpy(signal.CopyTree(str(OS & 'tau1_matched && tau2_matched' & loose_taus)))
# same sign
data_selected = asrootpy(data.CopyTree(str(SS & loose_taus)))

factory.AddSignalTree(signal_selected, 1.)
factory.AddBackgroundTree(data_selected, 1.)

variables = [
    ('tau1_BDTJetScore','F', (0, 1)),
    ('tau2_BDTJetScore','F', (0, 1)),
    ('sphericity','F', (0, 1)),
    ('sphericity_boosted','F', (0, 1)),
    ('aplanarity','F', (0, .5)),
    ('aplanarity_boosted','F', (0, .5)),
    ('tau1_x', 'F', (0, 1)),
    ('tau2_x', 'F', (0, 1)),
    ('cos_theta_tau1_tau2', 'F', (-1, 1)),
    ('jet1_jvtxf', 'F', (-1, 1)),
    ('jet2_jvtxf', 'F', (-1, 1)),
    ('dEta_jets', 'F', (0.1, 6)),
    ('dEta_jets_boosted', 'F', (0.1, 6))
]


for var, type, range in variables:

    canvas.Clear()
    sig = Hist(100, *range)
    bkg = Hist(100, *range)
    signal_selected.Draw(var, hist=sig)
    data_selected.Draw(var, hist=bkg)
    sig /= sum(sig)
    bkg /= sum(bkg)
    sig.format = 'hist'
    bkg.SetXTitle(var)
    bkg.SetMinimum(0)
    bkg.Draw()
    sig.Draw('same')
    canvas.SaveAs('plots/%s.eps' % var)
    canvas.SaveAs('plots/%s.png' % var)
    factory.AddVariable(var, type)

factory.PrepareTrainingAndTestTree(ROOT.TCut(''), "NormMode=EqualNumEvents")

params = ':'.join([
    'UseYesNoLeaf=False',
    'NTrees=10',
    'MaxDepth=100000',
    'BoostType=AdaBoost',
    'AdaBoostBeta=0.2',
    'SeparationType=GiniIndex',
    'PruneMethod=NoPruning',
    'PruneStrength=-1',
    'PruneBeforeBoost=False',
    'PruningValFraction=0.50',
    'NNodesMax=100000',
    'nCuts=20'
])

factory.BookMethod(ROOT.TMVA.Types.kBDT, "BDT", str(":".join(["!H","V",params])))

factory.TrainAllMethods()
factory.TestAllMethods()
factory.EvaluateAllMethods()
output.Close()
data_file.Close()
signal_file.Close()
tmp.Close()
os.unlink('tmp.root')

import ROOT
from atlastools.units import GeV
from atlastools import datasets
import os
from math import sqrt


ROOT.gSystem.Load(os.path.join(os.path.dirname(__file__), "libMissingMassCalculator.so"))
MMC = ROOT.MissingMassCalculator()
MMC.SetAlgorithmVersion(1)
MMC.SetNiterFit1(20)
MMC.SetNiterFit2(30)
#MMC.SetNiterFit3(10)
MMC.SetNsigmaMETscan(3.0)
MMC.SetApplyMassScale(0)


"""
This module contains utility functions for using libMissingMassCalculator.so
"""

def mass(tau1, tau2,
         METx, METy, sumET,
         tau2_lep_type = -1,
         method = 1):
    """
    Missing mass calculation
    returns the most likely mass
    """
    vis_tau1 = ROOT.TLorentzVector()
    # 1 prong
    if tau1.numTrack <= 1:
        tau1_decay_type = 10
        vis_tau1.SetPtEtaPhiM(tau1.pt/GeV, tau1.eta, tau1.phi, 0.8)
    # 3 prongs
    else:
        tau1_decay_type = 30
        vis_tau1.SetPtEtaPhiM(tau1.pt/GeV, tau1.eta, tau1.phi, 1.2)

    vis_tau2 = ROOT.TLorentzVector()

    if tau2_lep_type == -1:
        # 1 prong
        if tau2.numTrack <= 1:
            tau2_decay_type = 10
            vis_tau2.SetPtEtaPhiM(tau2.pt/GeV, tau2.eta, tau2.phi, 0.8)
        # 3 prongs
        else:
            tau2_decay_type = 30
            vis_tau2.SetPtEtaPhiM(tau2.pt/GeV, tau2.eta, tau2.phi, 1.2)

    # electron
    elif tau2_lep_type == 0:
        tau2_decay_type = 0
        vis_tau2.SetPtEtaPhiM(tau2.pt/GeV, tau2.eta, tau2.phi, 0.000511)

    # muon
    elif tau2_lep_type == 1:
        tau2_decay_type = 0
        vis_tau2.SetPtEtaPhiM(tau2.pt/GeV, tau2.eta, tau2.phi, 0.105658)

    else:
        raise ValueError('tau2_lep_type in missingmass.mass() should be -1 (had), 0 (electron) or 1 (muon). It is ' + str(tau2_lep_type))


    MMC.SetVisTauVec(0, vis_tau1)
    MMC.SetVisTauVec(1, vis_tau2)
    MMC.SetVisTauType(0, tau1_decay_type)
    MMC.SetVisTauType(1, tau2_decay_type)

    """
    jetvec = ROOT.vector("TLorentzVector")()
    jetP4 = ROOT.TLorentzVector()
    jet_SumEt = 0.

    for jet in jets:
        jetP4.SetPtEtaPhiM(jet.pt/GeV, jet.eta, jet.phi, jet.m/GeV);
        jet_SumEt += jetP4.Et();
        jetvec.push_back(jetP4);
    """

    """
    MMC_SumEt = sumET/GeV
    MMC_SumEt -= tau1.Et/GeV
    MMC_SumEt -= tau2.Et/GeV
    MMC_SumEt -= jet_SumEt
    """

    met_vec = ROOT.TVector2(METx/GeV, METy/GeV)
    MMC.SetMetVec(met_vec)

    MET_res = 6.14 + 0.5 * sqrt(abs(sumET) / GeV) # sumET can be negative!!
    MMC.SetMetScanParams(0.0, MET_res, MET_res)

    """
    if len(jets) > 0:
       MMC.SetMetScanParamsJets(jetvec)
    MMC.SetMetScanParamsUE(MMC_SumEt,  0.0, int(datatype == datasets.MC)) # data_type int : should be 1 for MC, 0 for data
    """

    MMC.RunMissingMassCalculator()
    MMC_mass = -1
    if MMC.GetFitStatus() == 1: # MMC output: 1=found solution; 0= no slution
        MMC_mass = MMC.GetFittedMass(method) # use 2 instead of 1 to remove spikes in output
    return MMC_mass

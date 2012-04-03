import ROOT
from atlastools.units import *
from rootpy.utils.classfactory import generate
from atlastools import datasets

ROOT.gSystem.Load("libMissingMassCalculator.so")
generate("vector<TLorentzVector>", "<vector>;TLorentzVector.h")

def mass(taus, jets, METx, METy, sumET, datatype):
    """
    Missing mass calculation
    returns the most likely mass
    """    
    jetvec = ROOT.vector("TLorentzVector")()
    jetP4 = ROOT.TLorentzVector()
    jet_SumEt = 0.

    for jet in jets:
        jetP4.SetPtEtaPhiM(jet.pt/GeV, jet.eta, jet.phi, jet.m/GeV);
        jet_SumEt += jetP4.Et();
        jetvec.push_back(jetP4);

    MMC_SumEt = sumET/GeV
    MMC_SumEt -= taus[0].Et/GeV
    MMC_SumEt -= taus[1].Et/GeV
    MMC_SumEt -= jet_SumEt
    
    VisTau0 = ROOT.TLorentzVector()
    if taus[0].seedCalo_numTrack <= 1:
        tau0_decay_type = 10
        VisTau0.SetPtEtaPhiM(taus[0].pt/GeV, taus[0].seedCalo_eta, taus[0].seedCalo_phi, 0.8)
    else:
        tau0_decay_type = 30
        VisTau0.SetPtEtaPhiM(taus[0].pt/GeV, taus[0].seedCalo_eta, taus[0].seedCalo_phi, 1.2)
    
    VisTau1 = ROOT.TLorentzVector()
    if taus[1].seedCalo_numTrack <= 1:
        tau1_decay_type = 10
        VisTau1.SetPtEtaPhiM(taus[1].pt/GeV, taus[1].seedCalo_eta, taus[1].seedCalo_phi, 0.8)
    else:
        tau1_decay_type = 30
        VisTau1.SetPtEtaPhiM(taus[1].pt/GeV, taus[1].seedCalo_eta, taus[1].seedCalo_phi, 1.2)
    
    met_vec = ROOT.TVector2(METx/GeV, METy/GeV)
    
    mmc = ROOT.MissingMassCalculator()
    mmc.SetMetVec(met_vec)
    mmc.SetVisTauVec(0, VisTau0)
    mmc.SetVisTauVec(1, VisTau1)
    mmc.SetVisTauType(0, tau0_decay_type)
    mmc.SetVisTauType(1, tau1_decay_type)
    mmc.SetMetScanParamsUE(MMC_SumEt,  0., int(datatype == datasets.MC)) # data_type int : should be 1 for MC, 0 for data
    misMassTest = mmc.RunMissingMassCalculator()
    output_fitstatus = mmc.GetFitStatus() # MMC output: 1=found solution; 0= no slution
    MMC_mass = 0.
    if output_fitstatus == 1:
        MMC_mass = mmc.GetFittedMass(2) # use 2 instead of 1 to remove spikes in output
    return MMC_mass

import ROOT
from atlastools.units import GeV
from atlastools import datasets
import os
from math import sqrt


class MMC(object):

    INITED = False

    def __init__(self, year, channel):

        if MMC.INITED:
            raise RuntimeError("do not create more than one MMC")
        MMC.INITED = True

        assert channel in ('lh', 'hh')
        assert year in (2011, 2012)

        if channel == 'hh' and year == 2012:
            from externaltools.bundle_2012 import MissingMassCalculator
            # tag 7
        else:
            from externaltools.bundle_2011 import MissingMassCalculator
            # tag 9
            # mmc tags are a mess...

        self.tool = ROOT.MissingMassCalculator()
        self.tool.SetAlgorithmVersion(1)
        self.tool.SetNiterFit1(20)
        self.tool.SetNiterFit2(30)
        #self.tool.SetNiterFit3(10)
        self.tool.SetNsigmaMETscan(3.0)
        self.tool.SetApplyMassScale(0)

        if channel == 'hh' and year == 2012:
            # temporary hack to reduce failure rate
            print "using larger MMC param space search..."
            self.tool.SetNiterFit2(40)
            self.tool.SetNsigmaMETscan(4.0)

    def mass(self,
            tau1, tau2,
            METx, METy, sumET,
            tau2_lep_type=-1,
            method=1,
            year=None,
            njets25=0):
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
            tau2_decay_type = 1
            vis_tau2.SetPtEtaPhiM(tau2.pt/GeV, tau2.eta, tau2.phi,0.10565836668)

        # muon
        elif tau2_lep_type == 1:
            tau2_decay_type = 0
            vis_tau2.SetPtEtaPhiM(tau2.pt/GeV, tau2.eta, tau2.phi, 0.000510999)

        else:
            raise ValueError(
                    'tau2_lep_type in missingmass.mass() should be '
                    '-1 (had), 0 (electron) or 1 (muon). It is ' + str(tau2_lep_type))

        self.tool.SetVisTauVec(0, vis_tau1)
        self.tool.SetVisTauVec(1, vis_tau2)
        self.tool.SetVisTauType(0, tau1_decay_type)
        self.tool.SetVisTauType(1, tau2_decay_type)

        if tau2_lep_type > -1:
            self.tool.SetSumEt(sumET)
            self.tool.SetNjet25(njets25)

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
        self.tool.SetMetVec(met_vec)

        if tau2_lep_type == -1:
            MET_res = 6.14 + 0.5 * sqrt(abs(sumET) / GeV) # sumET can be negative!!
            self.tool.SetMetScanParams(0.0, MET_res, MET_res)

        """
        if len(jets) > 0:
           MMC.SetMetScanParamsJets(jetvec)
        MMC.SetMetScanParamsUE(MMC_SumEt,  0.0, int(datatype == datasets.MC)) # data_type int : should be 1 for MC, 0 for data
        """

        self.tool.RunMissingMassCalculator()
        MMC_mass = -1
        MMC_resonance = ROOT.TLorentzVector(0, 0, 0, 0)
        MMC_met = ROOT.TVector2(0, 0)
        if self.tool.GetFitStatus() == 1: # MMC output: 1=found solution; 0= no slution
            MMC_mass = self.tool.GetFittedMass(method) # use 2 instead of 1 to remove spikes in output
            MMC_resonance = self.tool.GetResonanceVec(method)
            MMC_met = self.tool.GetFittedMetVec(method)
        return MMC_mass, MMC_resonance, MMC_met

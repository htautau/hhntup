from . import log; log = log[__name__]
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

        if channel == 'hh':
            log.info("setting up MMC for hadhad")
            from externaltools.bundle_2011 import MissingMassCalculator
            # tag 7 for hh for both years
        else:
            log.info("setting up MMC for lephad")
            from externaltools.bundle_2012 import MissingMassCalculator
            # tag 9 for lh for both years

        self.tool = ROOT.MissingMassCalculator('JERProviderPlots.root')

        if channel == 'hh':
            log.info("applying hadhad MMC settings")
            self.tool.SetAlgorithmVersion(1)
            self.tool.SetNiterFit2(40)
            self.tool.SetNsigmaMETscan(4.0)

    def mass(self,
             tau1, tau2,
             METx, METy, sumET,
             njets,
             tau2_lep_type=None,
             method=1):
        """
        Missing mass calculation
        returns the most likely mass
        """
        vis_tau1 = ROOT.TLorentzVector()
        # 1 prong
        if tau1.numTrack <= 1:
            tau1_decay_type = 10
            mvis = 0.8
        # 3 prongs
        else:
            tau1_decay_type = 30
            mvis = 1.2

        vis_tau1.SetPtEtaPhiM(
            tau1.fourvect.Pt() / GeV,
            tau1.fourvect.Eta(),
            tau1.fourvect.Phi(),
            mvis)

        vis_tau2 = ROOT.TLorentzVector()

        # hadronic tau
        if tau2_lep_type is None:
            # 1 prong
            if tau2.numTrack <= 1:
                tau2_decay_type = 10
                mvis = 0.8
            # 3 prongs
            else:
                tau2_decay_type = 30
                mvis = 1.2

            vis_tau2.SetPtEtaPhiM(
                tau2.fourvect.Pt() / GeV,
                tau2.fourvect.Eta(),
                tau2.fourvect.Phi(),
                mvis)

        # electron
        elif tau2_lep_type == 0:
            tau2_decay_type = 0
            vis_tau2.SetPtEtaPhiM(
                tau2.fourvect.Pt() / GeV,
                tau2.fourvect.Eta(),
                tau2.fourvect.Phi(),
                0.000510999)

        # muon
        elif tau2_lep_type == 1:
            tau2_decay_type = 1
            vis_tau2.SetPtEtaPhiM(
                tau2.fourvect.Pt() / GeV,
                tau2.fourvect.Eta(),
                tau2.fourvect.Phi(),
                0.105658367)

        else:
            raise ValueError(
                'tau2_lep_type in missingmass.mass() should be '
                '-1 (had), 0 (electron) or 1 (muon). It is %s' % tau2_lep_type)

        self.tool.SetVisTauVec(0, vis_tau1)
        self.tool.SetVisTauVec(1, vis_tau2)
        self.tool.SetVisTauType(0, tau1_decay_type)
        self.tool.SetVisTauType(1, tau2_decay_type)

        met_vec = ROOT.TVector2(METx / GeV, METy / GeV)
        self.tool.SetMetVec(met_vec)
        self.tool.SetSumEt(sumET / GeV)
        self.tool.SetNjet25(njets)

        MET_res = 6.14 + 0.5 * sqrt(abs(sumET / GeV)) # sumET can be negative!!
        self.tool.SetMetScanParams(0., MET_res, MET_res)

        self.tool.RunMissingMassCalculator()
        MMC_mass = -1
        MMC_resonance = ROOT.TLorentzVector(0, 0, 0, 0)
        MMC_met = ROOT.TVector2(0, 0)
        if self.tool.GetFitStatus() == 1:
            # MMC output: 1=found solution; 0= no solution
            MMC_mass = self.tool.GetFittedMass(method)
            # use method 2 instead of 1 to remove spikes in output
            MMC_resonance = self.tool.GetResonanceVec(method)
            MMC_met = self.tool.GetFittedMetVec(method)

        return MMC_mass, MMC_resonance, MMC_met

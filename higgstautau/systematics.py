"""
Adapted from the Example.C in MissingETUtility/macros
"""
# stdlib Python imports
import os
import math
from math import sin, sqrt, pow

# local imports
from . import tauid
from . import log; log = log[__name__]

from atlastools import utils
from atlastools import datasets

from rootpy.tree.filtering import EventFilter

# ROOT imports
import ROOT

# ATLAS tools imports
from externaltools import MissingETUtility
from externaltools import MuonMomentumCorrections
from externaltools import JetUncertainties
from externaltools import JetResolution
from externaltools import egammaAnalysisUtils
from externaltools import TauCorrUncert as TCU
from externaltools.bundle_2012 import JVFUncertaintyTool as JVFUncertaintyTool2012

# MissingETUtility
from ROOT import METUtility
from ROOT import METUtil
from ROOT import MissingETTags
# JetUncertainties
from ROOT import MultijetJESUncertaintyProvider
# JVF
from ROOT import JVFUncertaintyTool
from ROOT import TLorentzVector # needed for JVF
from rootpy import stl
VectorTLorentzVector = stl.vector("TLorentzVector")
# JetResolution
from ROOT import JERProvider
# egammaAnalysisUtils
from ROOT import eg2011
# MuonMomentumCorrections
from ROOT import MuonSmear
# tau corrections and uncertainties
from ROOT import TauCorrUncert


class ObjectSystematic(object):

    def __init__(self, sys_util):

        self.sys_util = sys_util
        self.verbose = sys_util.verbose
        self.year = sys_util.year
        self.channel = sys_util.channel


class JetSystematic(ObjectSystematic):

    def __init__(self, is_up, **kwargs):

        self.is_up = is_up # up or down variation
        super(JetSystematic, self).__init__(**kwargs)

    @staticmethod
    def set(f):

        def wrapper(self, event):

            if self.verbose:
                print "=" * 20
                print "JETS BEFORE:"
                for jet in event.jets:
                    print jet.pt, '\t',jet.jvtxf
                print "-" * 20

            for jet in event.jets:
                f(self, jet, event)

            if self.verbose:
                print "JETS AFTER:"
                for jet in event.jets:
                    print jet.pt, '\t',jet.jvtxf

        return wrapper


class JES(JetSystematic):
    # All 14 NPs in the 2013 recommendation
    all_nps = [
        'Statistical',
        'Modelling',
        'Detector',
        'Mixed',
        'EtaModelling',
        'EtaMethod',
        'PURho',
        'PUPt',
        'PUNPV',
        'PUMu',
        'FlavComp',
        'FlavResp',
        'BJet',
        'NonClosure',
    ]

    def __init__(self, is_up, np=None, **kwargs):

        # *** Set up the uncertainty tools ***
        # Tag assumed: JetUncertainties-00-05-09-02
        super(JES, self).__init__(is_up, **kwargs)

        assert np in self.all_nps and hasattr(self, '_'+np)
        self.np = np
        self.run_np = getattr(self, '_' + np)

        self.jes_tool = None
        if self.year == 2011:
            self.jes_tool = MultijetJESUncertaintyProvider(
                "JES_2011/Final/MultijetJES_2011.config",
                "JES_2011/Final/InsituJES2011_17NP_ByCategory.config",
                "AntiKt4TopoLC",
                "MC11c",
                JetUncertainties.RESOURCE_PATH)
        elif self.year == 2012:
            self.jes_tool = MultijetJESUncertaintyProvider(
                "JES_2012/Moriond2013/MultijetJES_2012.config",
                #"JES_2012/Moriond2013/InsituJES2012_AllNuisanceParameters.config",
                "JES_2012/Moriond2013/InsituJES2012_20NP_ByCategory.config",
                "AntiKt4TopoLC",
                "MC12a",
                JetUncertainties.RESOURCE_PATH)
        else:
            raise ValueError('No JES defined for year %d' % year)

    def _Statistical(self, jet, event):
        # HACK: use only stat1 for now
        uncs = [0.]
        #for i in range(3):
        uncs[0] = self.jes_tool.getRelUncertComponent(
            "EffectiveNP_Statistical1", jet.pt, jet.eta)
        return uncs

    def _Modelling(self, jet, event):
        # HACK: use only modelling1 for now
        uncs = [0.]
        #for i in range(4):
        uncs[0] = self.jes_tool.getRelUncertComponent(
            "EffectiveNP_Modelling1", jet.pt, jet.eta)
        return uncs

    def _Detector(self, jet, event):
        # HACK: use only detector1 for now
        uncs = [0.]
        #for i in range(3):
        uncs[0] = self.jes_tool.getRelUncertComponent(
            "EffectiveNP_Detector1", jet.pt, jet.eta)
        return uncs

    def _Mixed(self, jet, event):
        uncs = [0., 0.]
        for i in range(2):
            uncs[i] = self.jes_tool.getRelUncertComponent("EffectiveNP_Mixed%d" % (i+1), jet.pt, jet.eta)
        return uncs

    def _EtaModelling(self, jet, event):
        return [ self.jes_tool.getRelUncertComponent("EtaIntercalibration_Modelling", jet.pt, jet.eta) ]

    def _EtaMethod(self, jet, event):
        return [ self.jes_tool.getRelUncertComponent("EtaIntercalibration_StatAndMethod", jet.pt, jet.eta) ]

    def _PURho(self, jet, event):
        return [ self.jes_tool.getRelPileupRhoTopology(jet.pt, jet.eta) ]

    def _PUPt(self, jet, event):
        return [ self.jes_tool.getRelPileupPtTerm(jet.pt, jet.eta, self.sys_util.nvtxjets, event.averageIntPerXing) ]

    def _PUNPV(self, jet, event):
        return [ self.jes_tool.getRelNPVOffsetTerm(jet.pt, jet.eta, self.sys_util.nvtxjets) ]

    def _PUMu(self, jet, event):
        return [ self.jes_tool.getRelMuOffsetTerm(jet.pt, jet.eta, event.averageIntPerXing) ]

    def _FlavComp(self, jet, event):
        uncs = [0.]
        is_b = jet.flavor_truth_dRminToB < 0.4
        if not is_b:
            uncs[0] = self.jes_tool.getRelFlavorCompUncert(jet.pt, jet.eta, True)
        return uncs

    def _FlavResp(self, jet, event):
        uncs = [0.]
        is_b = jet.flavor_truth_dRminToB < 0.4
        if not is_b:
            uncs[0] = self.jes_tool.getRelFlavorResponseUncert(jet.pt, jet.eta)
        return uncs

    def _BJet(self, jet, event):
        uncs = [0.]
        is_b = jet.flavor_truth_dRminToB < 0.4
        if not is_b:
            uncs[0] = self.jes_tool.getRelBJESUncert(jet.pt, jet.eta)
        return uncs

    def _NonClosure(self, jet, event):
        uncs = [0., 0.]
        uncs[0] = self.jes_tool.getRelUncertComponent("SingleParticle_HighPt", jet.pt, jet.eta)
        uncs[1] = self.jes_tool.getRelUncertComponent("RelativeNonClosure_Pythia8", jet.pt, jet.eta)
        # TODO: If AFII, use "RelativeNonClosure_AFII" instead
        return uncs

    @JetSystematic.set
    def run(self, jet, event):
        # Safest to assume nothing about the uncertainties on soft jets.
        # These will go into SoftJets anyhow, and so the JES systematics
        # aren't used.
        shift = 0.
        # OLD 2011 code
        #if self.year == 2011:
        #     if jet.pt > 20e3 and jet.pt < 7000e3 and abs(jet.eta) < 4.5:
        #
        #         # delta R cut needed to apply close-by jets uncertainty
        #         drmin=9999
        #         for otherjet in event.jets:
        #             if otherjet.emscale_pt > 7000:
        #                 if jet.index != otherjet.index:
        #                     dr = utils.dR(jet.eta, jet.phi,
        #                                   otherjet.eta,
        #                                   otherjet.phi)
        #                     if dr < drmin:
        #                         drmin = dr
        #
        #         # TODO: shift is symmetric (is_up argument is not needed)
        #         if self.is_up:
        #             shift = self.jes_tool.getRelUncert(
        #                     jet.pt,
        #                     jet.eta,
        #                     drmin,
        #                     True, # is up
        #                     self.sys_util.nvtxjets,
        #                     event.averageIntPerXing,
        #                     False) # is b jet
        #         else:
        #             shift = -1 * self.jes_tool.getRelUncert(
        #                     jet.pt,
        #                     jet.eta,
        #                     drmin,
        #                     False, # is up
        #                     self.sys_util.nvtxjets,
        #                     event.averageIntPerXing,
        #                     False) # is b jet
        if jet.pt > 15e3 and jet.pt < 7000e3 and abs(jet.eta) < 4.5:

            uncs = self.run_np(jet, event)
            if len(uncs) > 1:
                unc = 0.
                for u in uncs:
                    unc += u**2.
                unc = unc**0.5
            else:
                unc = uncs[0]

            if abs(max(uncs)) < abs(min(uncs)):
                unc = -unc; # set total to sign of max
            if self.is_up:
                shift = unc;
            else:
                shift = -unc;
        jet.pt *= 1. + shift


class JVF(JetSystematic):

    def __init__(self, is_up, JVFcutNominal=0.5, **kwargs):
        self.JVFcutNominal = JVFcutNominal

        # Tag assumed: JVFUncertaintyTool-00-00-03
        self.jvf_tool = JVFUncertaintyTool("AntiKt4LCTopo")
        super(JVF, self).__init__(is_up, **kwargs)

    @JetSystematic.set
    def run(self, jet, event):
        """ From Melbourne's framework...:
        TLorentzVector jet;
        TLorentzVector aux_trueJet;
        std::vector<TLorentzVector> trueJets;

        for (int k=0; k < ao->truejet.n; k++)
        {
            aux_trueJet.SetPtEtaPhiM(   ao->truejet.pt->at(k),
                                        ao->truejet.eta->at(k),
                                        ao->truejet.phi->at(k),
                                        ao->truejet.m->at(k));
            if (aux_trueJet.Pt()<= 10000) continue;
            trueJets.push_back(aux_trueJet);
        }


        for (int j=0; j < ao->jet.n; j++)
        {
            jet.SetPtEtaPhiM(ao->jet.pt->at(j),  ao->jet.eta->at(j),  ao->jet.phi->at(j),  ao->jet.m->at(j));

            if (jet.Pt()>= 50000. || fabs(jet.Eta()) >= 2.4) continue;

            //Verify is the jet is classified as a PU or HS jet
            bool isPU = jvf_uncertainty_tool->isPileUpJet(jet, trueJets);

            double eta_det = ao->jet.constscale_eta->at(j);
            double JVFcutNominal = 0.5;
            float jvf_cut_sys = jvf_uncertainty_tool->getJVFcut(JVFcutNominal, isPU, jet.Pt(), eta_det, flag);

            //change variation from the jvf_cut to the jet jvf value
            //easier to implement in the current framework, where the jvf_cut is always fixed to 0.5
            float jvf_cut_diff = jvf_cut_sys - JVFcutNominal;
            ao->jet.jvtxf->at(j)  =  ao->jet.jvtxf->at(j) - jvf_cut_diff;
        }
        """
        # JVF is only used in a certain range, so only correct for those
        if jet.pt < 50e3 and abs(jet.eta) < 2.4:
            truejets = VectorTLorentzVector()
            truejets_cache = []
            for truejet in event.truejets:
                if truejet.pt > 10e3:
                    t = TLorentzVector()
                    t.SetPtEtaPhiM(truejet.pt, truejet.eta, truejet.phi, truejet.m)
                    truejets.push_back(t)
                    truejets_cache.append(t)

            j = TLorentzVector()
            j.SetPtEtaPhiM(jet.pt, jet.eta, jet.phi, jet.m)
            isPU = self.jvf_tool.isPileUpJet(j, truejets)
            jvf_cut_sys = self.jvf_tool.getJVFcut(self.JVFcutNominal, isPU, jet.pt, jet.constscale_eta, self.is_up)
            jvf_cut_diff = jvf_cut_sys - self.JVFcutNominal
            jet.jvtxf -= jvf_cut_diff

class JER(JetSystematic):

    def __init__(self, is_up, **kwargs):

        # Tag assumed: JetResolution-01-00-00
        if self.year == 2011:
            self.jer_tool = JERProvider(
                "AntiKt4LCTopoJES", "Truth",
                JetResolution.get_resource('JERProviderPlots_2011.root'))
            self.jer_tool.is7TeV(True)
        else:
            self.jer_tool = JERProvider(
                "AntiKt4LCTopoJES", "Truth",
                JetResolution.get_resource('JERProviderPlots.root'))
        self.jer_tool.init()
        # Note on use of ROOT random number generators:
        # TRandom and TRandom2 have many documented deficiencies.
        # TRandom3 is generally considered safely usable.
        # Also note that ROOT's gRandom calls TRandom3.
        self.jetrandom = ROOT.TRandom3()
        super(JER, self).__init__(is_up, **kwargs)

    @JetSystematic.set
    def run(self, jet, event):
        """
        Note: The JERDown shift is essentially meaningless.
        If one is smearing central values, then there is an alternate
        definition, i.e. from r16:

        S = self.jerTool.getRelResolutionData(pt/1e3,eta)
        SMC = self.jerTool.getRelResolutionMC(pt/1e3,eta)
        U = self.jerTool.getResolutionUncert(pt/1e3,eta)
        smearingFactorMC = sqrt( S*S - SMC*SMC )
        smearingFactorSystUp = sqrt( (S+U)*(S+U) - SMC*SMC )
        smearingFactorSystDown = (S-U > SMC) ? sqrt( (S+U)*(S+U) - SMC*SMC ) : 0

        float jerShift = jetRandom.Gaus(1,smearingFactorMC)
        float jerShiftUp = jetRandom.Gaus(1,smearingFactorSystUp)/jerShift
        float jerShiftDown = jetRandom.Gaus(1,smearingFactorSystDown)/jerShift

        jet_smeared_pt = pt*jerShift
        jerUp.push_back(jerShiftUp-1)
        jerDown.push_back(jerShiftDown-1)

        This means that we smear the MC jets to match the resolution in data
        for central values, or the resolution +/- uncertainty.
        The standard practice is only to use res + uncertainty.
        """
        shift = 0.
        # Allowable range is > 10 GeV, but anything below 20 enters SoftJets
        if jet.pt > 20e3 and jet.pt < 5000e3:
            pt = jet.pt
            eta = jet.eta
            if abs(eta) > 4.5:
                eta = 4.49 if eta > 0 else -4.49

            S = self.jer_tool.getRelResolutionMC(pt/1e3,eta)
            U = self.jer_tool.getResolutionUncert(pt/1e3,eta)
            smearingFactorSyst = sqrt(pow(S+U,2)-pow(S,2))

            # You can set the seed however you like, but if reproducibility
            # is required, setting it to something like object phi ensures
            # a good mix of randomness and reproducibility.
            self.jetrandom.SetSeed(int(1.e5 * abs(jet.phi)))
            shift = self.jetrandom.Gaus(0, smearingFactorSyst)

        if not self.is_up:
            shift *= -1
        jet.pt *= 1. + shift


class TauIDSystematic(ObjectSystematic):

    def __init__(self, is_up, **kwargs):

        self.is_up = is_up # up or down variation
        super(TauIDSystematic, self).__init__(**kwargs)

    @staticmethod
    def set(f):

        def wrapper(self, event):

            if self.verbose:
                print "=" * 20
                print "TAUS BEFORE:"
                for tau in event.taus:
                    print "score: %.4f loose: %d medium: %d tight: %d" % (
                        tau.BDTJetScore,
                        tau.JetBDTSigLoose,
                        tau.JetBDTSigMedium,
                        tau.JetBDTSigTight)
                print "-" * 20

            for tau in event.taus:
                f(self, tau, event)

            if self.verbose:
                print "TAUS AFTER:"
                for tau in event.taus:
                    print "score: %.4f loose: %d medium: %d tight: %d" % (
                        tau.BDTJetScore,
                        tau.JetBDTSigLoose,
                        tau.JetBDTSigMedium,
                        tau.JetBDTSigTight)

        return wrapper


class TauBDT(TauIDSystematic):
    """
    Currently only valid for 2011 MC
    """
    @TauIDSystematic.set
    def run(self, tau, event):

        high, low = tauid.uncertainty(
            tau.BDTJetScore, tau.pt, tau.numTrack,
            event.number_of_good_vertices)
        if self.is_up:
            tau.BDTJetScore = high
        else:
            tau.BDTJetScore = low


class TauSystematic(ObjectSystematic):

    def __init__(self, is_up, **kwargs):

        self.is_up = is_up # up or down variation
        super(TauSystematic, self).__init__(**kwargs)

    @staticmethod
    def set(f):

        def wrapper(self, event):

            if self.verbose:
                print "=" * 20
                print "TAUS BEFORE:"
                for tau in event.taus:
                    print tau.pt
                print "-" * 20

            for tau in event.taus:
                f(self, tau, event)

            if self.verbose:
                print "TAUS AFTER:"
                for tau in event.taus:
                    print tau.pt

        return wrapper


class TES(TauSystematic):

    def __init__(self, is_up,
                 np=TauCorrUncert.TESUncertainty.FINAL,
                 infile='TES/mc12_p1344_medium.root',
                 datatype=None,
                 matched_state=None,
                 **kwargs):

        super(TES, self).__init__(is_up, **kwargs)

        self.np = np
        self.datatype = datatype
        self.matched_state = matched_state
        if self.year == 2011:
            from externaltools.bundle_2011 import TESUncertaintyProvider as TESP
            from ROOT import TESUncertaintyProvider
            self.tes_tool = TESUncertaintyProvider(
                    os.path.normpath(TESP.RESOURCE_PATH), '', 'mc11')
        elif self.year == 2012:
            # TODO use medium and tight?
            self.tes_tool = TauCorrUncert.TESUncertainty(
                    TCU.get_resource(infile))
        else:
            raise ValueError('No TES defined for year %d' % self.year)

    @TauSystematic.set
    def run(self, tau, event):

        # need nominal pT for trigger efficiency correction
        tau._pt_nominal = tau.pt

        if self.matched_state is True and not tau.matched:
            return
        elif self.matched_state is False and tau.matched:
            return

        pt = tau.pt
        eta = tau.eta
        nProng = tau.nProng
        # TODO: include 2011 TES in TauCorrUncert and use MeV there also!!!
        t = TauCorrUncert.TESUncertainty
        if self.np == t.OTHERS and self.datatype in (datasets.DATA, datasets.EMBED):
            shift = 0.
            for np in [t.SHOWERMODEL, t.UE, t.DM]: # exclude CLOSURE manually
                shift += self.tes_tool.GetTESUncertainty(pt, eta, nProng, np)**2.
            shift = sqrt(shift)
        else:
            shift = self.tes_tool.GetTESUncertainty(pt, eta, nProng, self.np)
        if shift < 0:
            shift = 0
        if not self.is_up:
            shift *= -1
        tau.pt *= 1. + shift


class ElectronSystematic(ObjectSystematic):

    def __init__(self, is_up, **kwargs):

        self.is_up = is_up # up or down variation
        super(ElectronSystematic, self).__init__(**kwargs)

    @staticmethod
    def set(f):

        def wrapper(self, event):

            if self.verbose:
                print "=" * 20
                print "ELECTRONS BEFORE:"
                for el in event.electrons:
                    print el.pt
                print "-" * 20

            for el in event.electrons:
                f(self, el, event)

            if self.verbose:
                print "ELECTRONS AFTER:"
                for el in event.electrons:
                    print el.pt

        return wrapper


class EES(ElectronSystematic):

    def __init__(self, is_up, datatype, **kwargs):

        # Tag assumed: egammaAnalysisUtils-00-02-76
        self.egamma_tool = eg2011.EnergyRescaler()
        self.egamma_tool.useDefaultCalibConstants("2011")
        self.datatype = datatype
        super(EES, self).__init__(is_up, **kwargs)

    @ElectronSystematic.set
    def run(self, el, event):

        # Correct the measured energies in data, and scale by systematic variations
        correction = 1.
        if self.datatype == datasets.DATA:
            correction = self.egammaTool.applyEnergyCorrectionMeV(
                    el.cl_eta,
                    el.cl_phi,
                    el.E,
                    el.cl_pt,
                    0,"ELECTRON") / el.E
        if self.is_up:
            shift = self.egammaTool.applyEnergyCorrectionMeV(
                el.cl_eta,
                el.cl_phi,
                el.E,
                el.cl_pt,
                2,"ELECTRON") / (correction * el.E) - 1
        else:
            shift = self.egammaTool.applyEnergyCorrectionMeV(
                el.cl_eta,
                el.cl_phi,
                el.E,
                el.cl_pt,
                1,"ELECTRON") / (correction * el.E) - 1
        el.pt *= 1. + shift


class EER(ElectronSystematic):

    def __init__(self, is_up, **kwargs):

        # Tag assumed: egammaAnalysisUtils-00-02-76
        self.egamma_tool = eg2011.EnergyRescaler()
        self.egamma_tool.useDefaultCalibConstants("2011")
        super(EER, self).__init__(is_up, **kwargs)

    @ElectronSystematic.set
    def run(self, el, event):

        self.egamma_tool.SetRandomSeed(int(1e5*abs(el.phi)))
        # Smear to match the data resolution, or by systematic variations
        smear_nominal = self.egamma_tool.getSmearingCorrectionMeV(el.cl_eta, el.E, 0, True)
        if self.is_up:
            smear = self.egamma_tool.getSmearingCorrectionMeV(el.cl_eta, el.E, 2, True)
        else:
            smear = self.egamma_tool.getSmearingCorrectionMeV(el.cl_eta, el.E, 1, True)
        el.pt *= 1. + ((smear - smear_nominal) / smear_nominal)


class Systematics(EventFilter):

    # default
    NONE = METUtil.None

    # electrons
    EES_UP = METUtil.EESUp
    EES_DOWN = METUtil.EESDown
    EER_UP = METUtil.EERUp
    EER_DOWN = METUtil.EERDown
    ELECTRON_TERMS = set([EES_UP, EES_DOWN, EER_UP, EER_DOWN])

    # photons
    PES_UP = METUtil.PESUp
    PES_DOWN = METUtil.PESDown
    PER_UP = METUtil.PERUp
    PER_DOWN = METUtil.PERDown
    PHOTON_TERMS = set([PES_UP, PES_DOWN, PER_UP, PER_DOWN])

    # taus
    TES_UP = METUtil.TESUp
    TES_DOWN = METUtil.TESDown

    TES_TRUE_UP = -30000
    TES_TRUE_DOWN = -30001
    TES_FAKE_UP = -30010
    TES_FAKE_DOWN = -30011

    TES_EOP_UP = -20000
    TES_EOP_DOWN = -20001
    TES_CTB_UP = -20010
    TES_CTB_DOWN = -20011
    TES_Bias_UP = -20020
    TES_Bias_DOWN = -20021
    TES_EM_UP = -20030
    TES_EM_DOWN = -20031
    TES_LCW_UP = -20040
    TES_LCW_DOWN = -20041
    TES_PU_UP = -20050
    TES_PU_DOWN = -20051
    TES_OTHERS_UP = -20060
    TES_OTHERS_DOWN = -20061

    TER_UP = METUtil.TERUp
    TER_DOWN = METUtil.TERDown

    TES_TERMS = set([
        TES_UP, TES_DOWN,
        TES_TRUE_UP, TES_TRUE_DOWN,
        TES_FAKE_UP, TES_FAKE_DOWN,
        TES_EOP_UP,
        TES_EOP_DOWN,
        TES_CTB_UP,
        TES_CTB_DOWN,
        TES_Bias_UP,
        TES_Bias_DOWN,
        TES_EM_UP,
        TES_EM_DOWN,
        TES_LCW_UP,
        TES_LCW_DOWN,
        TES_PU_UP,
        TES_PU_DOWN,
        TES_OTHERS_UP,
        TES_OTHERS_DOWN])

    TAU_TERMS = set([TER_UP, TER_DOWN]) | TES_TERMS

    TAUBDT_UP = -100
    TAUBDT_DOWN = -101

    # jets
    JES_UP = METUtil.JESUp
    JES_DOWN = METUtil.JESDown

    JES_Statistical_UP = -10000
    JES_Statistical_DOWN = -10001
    JES_Modelling_UP = -10010
    JES_Modelling_DOWN = -10011
    JES_Detector_UP = -10020
    JES_Detector_DOWN = -10021
    JES_Mixed_UP = -10030
    JES_Mixed_DOWN = -10031
    JES_EtaModelling_UP = -10040
    JES_EtaModelling_DOWN = -10041
    JES_EtaMethod_UP = -10050
    JES_EtaMethod_DOWN = -10051
    JES_PURho_UP = -10060
    JES_PURho_DOWN = -10061
    JES_PUPt_UP = -10070
    JES_PUPt_DOWN = -10071
    JES_PUNPV_UP = -10080
    JES_PUNPV_DOWN = -10081
    JES_PUMu_UP = -10090
    JES_PUMu_DOWN = -10091
    JES_FlavComp_UP = -10100
    JES_FlavComp_DOWN = -10101
    JES_FlavResp_UP = -10110
    JES_FlavResp_DOWN = -10111
    JES_BJet_UP = -10120
    JES_BJet_DOWN = -10121
    JES_NonClosure_UP = -10130
    JES_NonClosure_DOWN = -10131

    JER_UP = METUtil.JERUp
    JER_DOWN = METUtil.JERDown # NOT USED!

    JVF_UP = -11000
    JVF_DOWN = -11001

    JES_TERMS = set([
        JES_UP,                 JES_DOWN,
        JES_Statistical_UP,     JES_Statistical_DOWN,
        JES_Modelling_UP,       JES_Modelling_DOWN,
        JES_Detector_UP,        JES_Detector_DOWN,
        JES_Mixed_UP,           JES_Mixed_DOWN,
        JES_EtaModelling_UP,    JES_EtaModelling_DOWN,
        JES_EtaMethod_UP,       JES_EtaMethod_DOWN,
        JES_PURho_UP,           JES_PURho_DOWN,
        JES_PUPt_UP,            JES_PUPt_DOWN,
        JES_PUNPV_UP,           JES_PUNPV_DOWN,
        JES_PUMu_UP,            JES_PUMu_DOWN,
        JES_FlavComp_UP,        JES_FlavComp_DOWN,
        JES_FlavResp_UP,        JES_FlavResp_DOWN,
        JES_BJet_UP,            JES_BJet_DOWN,
        JES_NonClosure_UP,      JES_NonClosure_DOWN,
        JVF_UP,                 JVF_DOWN,
        JER_UP,                 JER_DOWN,
    ])

    # muons
    MERID_UP = METUtil.MERIDUp
    MERID_DOWN = METUtil.MERIDDown
    MERMS_UP = METUtil.MERMSUp
    MERMS_DOWN = METUtil.MERMSDown
    MES_UP = METUtil.MESUp
    MES_DOWN = METUtil.MESDown
    MUON_TERMS = set([
        MERID_UP, MERID_DOWN,
        MERMS_UP, MERMS_DOWN,
        MES_UP, MES_DOWN])

    # tracks
    TRKES_UP = METUtil.TrkESUp
    TRKES_DOWN = METUtil.TrkESDown
    TRKER_UP = METUtil.TrkERUp
    TRKER_DOWN = METUtil.TrkERDown

    # clusters
    CES_UP = METUtil.CESUp
    CES_DOWN = METUtil.CESDown
    CER_UP = METUtil.CERUp
    CER_DOWN = METUtil.CERDown

    # deprecated 1.1.0
    #ALLCLUSTERS_UP = METUtil.AllClustersUp
    #ALLCLUSTERS_DOWN = METUtil.AllClustersDown

    # soft terms
    MET_RESOSOFTTERMS_PTHARD_UP = METUtil.ResoSoftTermsUp_ptHard
    MET_RESOSOFTTERMS_PTHARD_DOWN = METUtil.ResoSoftTermsDown_ptHard

    MET_RESOSOFTTERMS_PTHARD_UPDOWN = METUtil.ResoSoftTermsUpDown_ptHard
    MET_RESOSOFTTERMS_PTHARD_DOWNUP = METUtil.ResoSoftTermsDownUp_ptHard

    MET_SCALESOFTTERMS_PTHARD_UP = METUtil.ScaleSoftTermsUp_ptHard
    MET_SCALESOFTTERMS_PTHARD_DOWN = METUtil.ScaleSoftTermsDown_ptHard

    MET_RESOSOFTTERMS_UP = METUtil.ResoSoftTermsUp
    MET_RESOSOFTTERMS_DOWN = METUtil.ResoSoftTermsDown

    MET_SCALESOFTTERMS_UP = METUtil.ScaleSoftTermsUp
    MET_SCALESOFTTERMS_DOWN = METUtil.ScaleSoftTermsDown
    MET_TERMS = set([
        MET_RESOSOFTTERMS_UP, MET_RESOSOFTTERMS_DOWN,
        MET_SCALESOFTTERMS_UP, MET_SCALESOFTTERMS_DOWN])

    # pileup (deprecated 1.1.0)
    #PILEUP_UP = METUtil.PileupUp
    #PILEUP_DOWN = METUtil.PileupDown

    def __init__(self,
            datatype,
            year,
            tree,
            channel='hh',
            terms=None,
            verbose=False,
            very_verbose=False,
            **kwargs):

        super(Systematics, self).__init__(**kwargs)

        self.systematics = []
        self.terms = set([])
        self.datatype = datatype
        self.year = year
        self.tree = tree
        self.channel = channel
        self.verbose = verbose
        self.very_verbose = very_verbose

        if terms is not None:
            # remove possible duplicates
            terms = set(terms)
            self.terms = terms
            for term in terms:
                systematic = None
                #if term == Systematics.JES_UP:
                #    systematic = JES(True, sys_util=self)
                #elif term == Systematics.JES_DOWN:
                #    systematic = JES(False, sys_util=self)
                if term == Systematics.JES_Statistical_UP:
                    systematic = JES(True, np="Statistical", sys_util=self)
                elif term == Systematics.JES_Statistical_DOWN:
                    systematic = JES(False,  np="Statistical", sys_util=self)
                elif term == Systematics.JES_Modelling_UP:
                    systematic = JES(True, np="Modelling", sys_util=self)
                elif term == Systematics.JES_Modelling_DOWN:
                    systematic = JES(False,  np="Modelling", sys_util=self)
                elif term == Systematics.JES_Detector_UP:
                    systematic = JES(True, np="Detector", sys_util=self)
                elif term == Systematics.JES_Detector_DOWN:
                    systematic = JES(False,  np="Detector", sys_util=self)
                elif term == Systematics.JES_Mixed_UP:
                    systematic = JES(True, np="Mixed", sys_util=self)
                elif term == Systematics.JES_Mixed_DOWN:
                    systematic = JES(False,  np="Mixed", sys_util=self)
                elif term == Systematics.JES_EtaModelling_UP:
                    systematic = JES(True, np="EtaModelling", sys_util=self)
                elif term == Systematics.JES_EtaModelling_DOWN:
                    systematic = JES(False,  np="EtaModelling", sys_util=self)
                elif term == Systematics.JES_EtaMethod_UP:
                    systematic = JES(True, np="EtaMethod", sys_util=self)
                elif term == Systematics.JES_EtaMethod_DOWN:
                    systematic = JES(False,  np="EtaMethod", sys_util=self)
                elif term == Systematics.JES_PURho_UP:
                    systematic = JES(True, np="PURho", sys_util=self)
                elif term == Systematics.JES_PURho_DOWN:
                    systematic = JES(False,  np="PURho", sys_util=self)
                elif term == Systematics.JES_PUPt_UP:
                    systematic = JES(True, np="PUPt", sys_util=self)
                elif term == Systematics.JES_PUPt_DOWN:
                    systematic = JES(False,  np="PUPt", sys_util=self)
                elif term == Systematics.JES_PUNPV_UP:
                    systematic = JES(True, np="PUNPV", sys_util=self)
                elif term == Systematics.JES_PUNPV_DOWN:
                    systematic = JES(False,  np="PUNPV", sys_util=self)
                elif term == Systematics.JES_PUMu_UP:
                    systematic = JES(True, np="PUMu", sys_util=self)
                elif term == Systematics.JES_PUMu_DOWN:
                    systematic = JES(False,  np="PUMu", sys_util=self)
                elif term == Systematics.JES_FlavComp_UP:
                    systematic = JES(True, np="FlavComp", sys_util=self)
                elif term == Systematics.JES_FlavComp_DOWN:
                    systematic = JES(False,  np="FlavComp", sys_util=self)
                elif term == Systematics.JES_FlavResp_UP:
                    systematic = JES(True, np="FlavResp", sys_util=self)
                elif term == Systematics.JES_FlavResp_DOWN:
                    systematic = JES(False,  np="FlavResp", sys_util=self)
                elif term == Systematics.JES_BJet_UP:
                    systematic = JES(True, np="BJet", sys_util=self)
                elif term == Systematics.JES_BJet_DOWN:
                    systematic = JES(False,  np="BJet", sys_util=self)
                elif term == Systematics.JES_NonClosure_UP:
                    systematic = JES(True, np="NonClosure", sys_util=self)
                elif term == Systematics.JES_NonClosure_DOWN:
                    systematic = JES(False,  np="NonClosure", sys_util=self)
                elif term == Systematics.JVF_UP:
                    systematic = JVF(True, sys_util=self)
                elif term == Systematics.JVF_DOWN:
                    systematic = JVF(False, sys_util=self)
                elif term == Systematics.JER_UP:
                    systematic = JER(True, sys_util=self)

                elif term == Systematics.TES_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.FINAL,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.FINAL,
                            infile='TES/mc12_p1344_medium_split.root')

                elif term == Systematics.TES_TRUE_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.FINAL,
                            infile='TES/mc12_p1344_medium_split.root',
                            matched_state=True)
                elif term == Systematics.TES_TRUE_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.FINAL,
                            infile='TES/mc12_p1344_medium_split.root',
                            matched_state=True)

                elif term == Systematics.TES_FAKE_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.FINAL,
                            infile='TES/mc12_p1344_medium_split.root',
                            matched_state=False)
                elif term == Systematics.TES_FAKE_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.FINAL,
                            infile='TES/mc12_p1344_medium_split.root',
                            matched_state=False)

                elif term == Systematics.TES_EOP_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.EOP,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_EOP_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.EOP,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_CTB_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.CTB,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_CTB_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.CTB,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_Bias_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.Bias,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_Bias_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.Bias,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_EM_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.EM,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_EM_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.EM,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_LCW_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.LCW,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_LCW_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.LCW,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_PU_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.PU,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_PU_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.PU,
                            infile='TES/mc12_p1344_medium_split.root')
                elif term == Systematics.TES_OTHERS_UP:
                    systematic = TES(True, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.OTHERS,
                            infile='TES/mc12_p1344_medium_split.root',
                            datatype=datatype)
                elif term == Systematics.TES_OTHERS_DOWN:
                    systematic = TES(False, sys_util=self,
                            np=TauCorrUncert.TESUncertainty.OTHERS,
                            infile='TES/mc12_p1344_medium_split.root',
                            datatype=datatype)

                elif term == Systematics.EES_UP:
                    systematic = EES(True, sys_util=self,
                        datatype=datatype)
                elif term == Systematics.EES_DOWN:
                    systematic = EES(False, sys_util=self,
                        datatype=datatype)
                elif term == Systematics.EER_UP:
                    systematic = EER(True, sys_util=self)
                elif term == Systematics.EER_DOWN:
                    systematic = EER(False, sys_util=self)
                elif term == Systematics.TAUBDT_UP:
                    systematic = TauBDT(True, sys_util=self)
                elif term == Systematics.TAUBDT_DOWN:
                    systematic = TauBDT(False, sys_util=self)
                elif term not in Systematics.MET_TERMS:
                    raise ValueError("systematic not supported")
                if systematic is not None:
                    self.systematics.append(systematic)

        # Initialise your METUtility object
        # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/MissingETUtility
        self.met_utility = METUtility()

        # configure
        if self.year == 2012:
            self.met_utility.configMissingET(
                year == 2012, # is 2012
                year == 2012) # is STVF

            # Turn on (off) the relevant MET terms
            # These are the terms required for MET_RefFinal(_BDTMedium)
            self.met_utility.defineMissingET(
                    True,  # RefEle
                    True,  # RefGamma
                    True,  # RefTau
                    True,  # RefJet
                    True,  # RefMuon
                    True,  # MuonTotal
                    True,  # Soft
                )

        # The threshold below which jets enter the SoftJets term (JES is not applied)
        #self.met_utility.setSoftJetCut(20e3)

        # Whether to use MUID muons (otherwise STACO).
        self.met_utility.setIsMuid(False)

        # Whether METUtility should scream at you over every little thing
        self.met_utility.setVerbosity(self.very_verbose)

        # Tag assumed: MuonMomentumCorrections-00-05-03
        # TODO: set the proper year here...
        #self.muonTool = MuonSmear.SmearingClass(
        #        "Data11","staco","pT","Rel17",
        #        MuonMomentumCorrections.RESOURCE_PATH)
        #self.muonTool.UseScale(1)
        #self.muonTool.UseImprovedCombine()

    def get_met(self):

        util = self.met_utility
        if self.year == 2011:
            if Systematics.MET_SCALESOFTTERMS_UP in self.terms:
                return util.getMissingET(METUtil.RefFinal, METUtil.ScaleSoftTermsUp)
            elif Systematics.MET_SCALESOFTTERMS_DOWN in self.terms:
                return util.getMissingET(METUtil.RefFinal, METUtil.ScaleSoftTermsDown)
            elif Systematics.MET_RESOSOFTTERMS_UP in self.terms:
                return util.getMissingET(METUtil.RefFinal, METUtil.ResoSoftTermsUp)
            elif Systematics.MET_RESOSOFTTERMS_DOWN in self.terms:
                return util.getMissingET(METUtil.RefFinal, METUtil.ResoSoftTermsDown)
            return util.getMissingET(METUtil.RefFinal)
        elif self.year == 2012:
            multisyst = METUtil.MultiSyst()
            if Systematics.MET_SCALESOFTTERMS_UP in self.terms:
                multisyst.setSyst(Systematics.MET_SCALESOFTTERMS_UP)
                return util.getMissingET(METUtil.RefFinal, multisyst)
            elif Systematics.MET_SCALESOFTTERMS_DOWN in self.terms:
                multisyst.setSyst(Systematics.MET_SCALESOFTTERMS_DOWN)
                return util.getMissingET(METUtil.RefFinal, multisyst)
            elif Systematics.MET_RESOSOFTTERMS_UP in self.terms:
                multisyst.setSyst(Systematics.MET_RESOSOFTTERMS_UP)
                return util.getMissingET(METUtil.RefFinal, multisyst)
            elif Systematics.MET_RESOSOFTTERMS_DOWN in self.terms:
                multisyst.setSyst(Systematics.MET_RESOSOFTTERMS_DOWN)
                return util.getMissingET(METUtil.RefFinal, multisyst)
            return util.getMissingET(METUtil.RefFinal)

    def passes(self, event):

        # reset the METUtility
        self.met_utility.reset()

        # Check for a good primary vertex
        # This is needed for jet and soft term systematics
        goodPV = False
        nvtxsoftmet = 0
        self.nvtxjets = 0

        # Most D3PDs contain the vx_type branch, but some don't.
        # Those which don't are most likely skimmed, and require at least 1
        # primary vertex for all events.
        # If your D3PD is skimmed in this way, then the goodPV (nTracks and z)
        # check should be applied to the first vertex in the collection.
        # Otherwise, you should ensure that the vx_type branch is available.
        for vertex in event.vertices:
            if vertex.type == 1 and vertex.nTracks > 2 and abs(vertex.z) < 200:
                goodPV = True
        if goodPV:
            for vertex in event.vertices:
                if vertex.nTracks > 2:
                    nvtxsoftmet += 1
                if vertex.nTracks > 1:
                    self.nvtxjets += 1

        # These set up the systematic "SoftTerms_ptHard"
        self.met_utility.setNvtx(nvtxsoftmet)

        # ResoSoftTerms uses gRandom for smearing.
        # Set the seed here however you like.
        if self.datatype in (datasets.DATA, datasets.EMBED):
            ROOT.gRandom.SetSeed(int(event.RunNumber * event.EventNumber))
        else:
            ROOT.gRandom.SetSeed(int(event.mc_channel_number * event.EventNumber))

        for systematic in self.systematics:
            systematic.run(event)

        if self.year == 2012:
            return self.passes_12(event)
        elif self.year == 2011:
            return self.passes_11(event)
        else:
            raise ValueError("No MET defined for year %d" % self.year)

    def passes_11(self, event):
        """
        JETS
        Always use setJetParameters since they may be recalibrated upstream
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/MissingETUtilityFAQ#If_I_recalibrate_correct_my_anal
        """
        if self.channel == 'hh':
            self.met_utility.setJetParameters(
                event.jet_pt,
                self.tree.jet_eta_original,
                self.tree.jet_phi_original,
                event.jet_E,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_wet,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_wpx,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_wpy,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_statusWord)

            #self.met_utility.setOriJetParameters(self.tree.jet_pt_original)

        if self.channel == 'lh':
            self.met_utility.setJetParameters(
                event.jet_pt,
                event.jet_eta,
                event.jet_phi,
                event.jet_PhiOriginEM,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_wet,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_wpx,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_wpy,
                event.jet_AntiKt4LCTopo_MET_BDTMedium_statusWord)

            #self.met_utility.setOriJetParameters(event.jet_EtaOriginEM)

        """ NEVER USE THIS since jets may be recalibrated upstream
        self.met_utility.setMETTerm(
            METUtil.RefJet,
            event.MET_RefJet_BDTMedium_etx,
            event.MET_RefJet_BDTMedium_ety,
            event.MET_RefJet_BDTMedium_sumet)
        """

        """
        ELECTRONS
        """
        if self.channel == 'hh':
            if self.terms & Systematics.ELECTRON_TERMS:
                self.met_utility.setElectronParameters(
                    event.el_pt,
                    event.el_eta,
                    event.el_phi,
                    event.el_MET_BDTMedium_wet,
                    event.el_MET_BDTMedium_wpx,
                    event.el_MET_BDTMedium_wpy,
                    event.el_MET_BDTMedium_statusWord)
            else:
                self.met_utility.setMETTerm(
                    METUtil.RefEle,
                    event.MET_RefEle_BDTMedium_etx,
                    event.MET_RefEle_BDTMedium_ety,
                    event.MET_RefEle_BDTMedium_sumet)

        if self.channel == 'lh':
            self.met_utility.setElectronParameters(
                event.el_cl_pt,
                event.el_eta,
                event.el_phi,
                event.el_MET_BDTMedium_wet,
                event.el_MET_BDTMedium_wpx,
                event.el_MET_BDTMedium_wpy,
                event.el_MET_BDTMedium_statusWord)


        if self.terms & Systematics.PHOTON_TERMS:
            self.met_utility.setPhotonParameters(
                event.ph_pt,
                event.ph_eta,
                event.ph_phi,
                event.ph_MET_BDTMedium_wet,
                event.ph_MET_BDTMedium_wpx,
                event.ph_MET_BDTMedium_wpy,
                event.ph_MET_BDTMedium_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefGamma,
                event.MET_RefGamma_BDTMedium_etx,
                event.MET_RefGamma_BDTMedium_ety,
                event.MET_RefGamma_BDTMedium_sumet)

        """
        MUONS
        """
        if self.channel == 'hh':
            if self.terms & Systematics.MUON_TERMS:
                self.met_utility.setMuonParameters(
                    event.mu_staco_pt, # or smeared pT
                    event.mu_staco_eta,
                    event.mu_staco_phi,
                    event.mu_staco_MET_BDTMedium_wet,
                    event.mu_staco_MET_BDTMedium_wpx,
                    event.mu_staco_MET_BDTMedium_wpy,
                    event.mu_staco_MET_BDTMedium_statusWord)

                # In this instance there is an overloaded version of
                # setExtraMuonParameters that accepts smeared pTs for spectro
                self.met_utility.setExtraMuonParameters(
                    event.mu_staco_ms_pt, # or smeared pT
                    event.mu_staco_ms_theta,
                    event.mu_staco_ms_phi)
            else:
                self.met_utility.setMETTerm(
                    METUtil.MuonTotal,
                    event.MET_Muon_Total_Staco_BDTMedium_etx,
                    event.MET_Muon_Total_Staco_BDTMedium_ety,
                    event.MET_Muon_Total_Staco_BDTMedium_sumet)

        if self.channel == 'lh':
            self.met_utility.setMuonParameters(
                event.mu_staco_etcone30, # or smeared pT
                event.mu_staco_eta,
                event.mu_staco_phi,
                event.mu_staco_MET_BDTMedium_wet,
                event.mu_staco_MET_BDTMedium_wpx,
                event.mu_staco_MET_BDTMedium_wpy,
                event.mu_staco_MET_BDTMedium_statusWord)

            # In this instance there is an overloaded version of
            # setExtraMuonParameters that accepts smeared pTs for spectro
            self.met_utility.setExtraMuonParameters(
                event.mu_staco_ptcone30, # or smeared pT
                event.mu_staco_ms_theta,
                event.mu_staco_ms_phi)


        # Note that RefMuon is not rebuilt from muons
        # -- it is a calorimeter term.
        self.met_utility.setMETTerm(
            METUtil.RefMuon,
            event.MET_RefMuon_Staco_BDTMedium_etx,
            event.MET_RefMuon_Staco_BDTMedium_ety,
            event.MET_RefMuon_Staco_BDTMedium_sumet)

        """
        TAUS
        """
        if self.channel == 'hh':
            if self.terms & Systematics.TAU_TERMS:
                self.met_utility.setTauParameters(
                    event.tau_pt,
                    event.tau_eta,
                    event.tau_phi,
                    event.tau_MET_BDTMedium_wet,
                    event.tau_MET_BDTMedium_wpx,
                    event.tau_MET_BDTMedium_wpy,
                    event.tau_MET_BDTMedium_statusWord)
            else:
                self.met_utility.setMETTerm(
                    METUtil.RefTau,
                    event.MET_RefTau_BDTMedium_etx,
                    event.MET_RefTau_BDTMedium_ety,
                    event.MET_RefTau_BDTMedium_sumet)

        if self.channel == 'lh':
            self.met_utility.setTauParameters(
                event.tau_pt,
                event.tau_eta,
                event.tau_phi,
                event.tau_MET_BDTMedium_wet,
                event.tau_MET_BDTMedium_wpx,
                event.tau_MET_BDTMedium_wpy,
                event.tau_MET_BDTMedium_statusWord)

        """
        self.met_utility.setMETTerm(
            METUtil.SoftJets,
            event.MET_SoftJets_BDTMedium_etx,
            event.MET_SoftJets_BDTMedium_ety,
            event.MET_SoftJets_BDTMedium_sumet)
        """

        self.met_utility.setMETTerm(
            METUtil.SoftTerms,
            event.MET_CellOut_BDTMedium_etx,
            event.MET_CellOut_BDTMedium_ety,
            event.MET_CellOut_BDTMedium_sumet)

        MET = self.get_met()

        if self.verbose:
            log.info("Recalculated MET: %.3f (original: %.3f)" % (
                     MET.et(), event.MET_RefFinal_BDTMedium_et))

        # update the MET with the shifted value
        self.tree.MET_etx_original = event.MET_RefFinal_BDTMedium_etx
        self.tree.MET_ety_original = event.MET_RefFinal_BDTMedium_ety
        self.tree.MET_et_original = event.MET_RefFinal_BDTMedium_et
        self.tree.MET_sumet_original = event.MET_RefFinal_BDTMedium_sumet
        self.tree.MET_phi_original = event.MET_RefFinal_BDTMedium_phi

        event.MET_RefFinal_BDTMedium_etx = MET.etx()
        event.MET_RefFinal_BDTMedium_ety = MET.ety()
        event.MET_RefFinal_BDTMedium_et = MET.et()
        event.MET_RefFinal_BDTMedium_sumet = MET.sumet()
        event.MET_RefFinal_BDTMedium_phi = MET.phi()

        return True

    def passes_12(self, event):

        # this must be put before setting the jets parameters
        self.met_utility.setJetPUcode(MissingETTags.JPU_JET_JVFCUT)

        """
        JETS
        Always use setJetParameters since they may be recalibrated upstream
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/MissingETUtilityFAQ#If_I_recalibrate_correct_my_anal
        """
        #if self.verbose:
        #    log.info(', '.join(map(str, self.tree.jet_phi_original)))
        #    log.info(', '.join(map(str, event.jet_phi)))
        self.met_utility.setJetParameters(
            event.jet_pt,
            self.tree.jet_eta_original,
            self.tree.jet_phi_original,
            event.jet_E,
            event.jet_AntiKt4LCTopo_MET_wet,
            event.jet_AntiKt4LCTopo_MET_wpx,
            event.jet_AntiKt4LCTopo_MET_wpy,
            event.jet_AntiKt4LCTopo_MET_statusWord)

        #self.met_utility.setOriJetParameters(event.jet_pt)

        """ NEVER USE THIS since jets may be recalibrated upstream
        self.met_utility.setMETTerm(
            METUtil.RefJet,
            event.MET_RefJet_STVF_etx,
            event.MET_RefJet_STVF_ety,
            event.MET_RefJet_STVF_sumet)
        """

        """
        ELECTRONS
        """
        if self.terms & Systematics.ELECTRON_TERMS:
            self.met_utility.setElectronParameters(
                event.el_pt,
                event.el_eta,
                event.el_phi,
                event.el_MET_wet,
                event.el_MET_wpx,
                event.el_MET_wpy,
                event.el_MET_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefEle,
                event.MET_RefEle_etx,
                event.MET_RefEle_ety,
                event.MET_RefEle_sumet)

        if self.terms & Systematics.PHOTON_TERMS:
            self.met_utility.setPhotonParameters(
                event.ph_pt,
                event.ph_eta,
                event.ph_phi,
                event.ph_MET_wet,
                event.ph_MET_wpx,
                event.ph_MET_wpy,
                event.ph_MET_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefGamma,
                event.MET_RefGamma_etx,
                event.MET_RefGamma_ety,
                event.MET_RefGamma_sumet)

        """
        MUONS
        """
        if self.terms & Systematics.MUON_TERMS:
            self.met_utility.setMuonParameters(
                event.mu_staco_pt, # or smeared pT
                event.mu_staco_eta,
                event.mu_staco_phi,
                event.mu_staco_MET_wet,
                event.mu_staco_MET_wpx,
                event.mu_staco_MET_wpy,
                event.mu_staco_MET_statusWord)

            # In this instance there is an overloaded version of
            # setExtraMuonParameters that accepts smeared pTs for spectro
            self.met_utility.setExtraMuonParameters(
                event.mu_staco_ms_pt, # or smeared pT
                event.mu_staco_ms_theta,
                event.mu_staco_ms_phi)
        else:
            self.met_utility.setMETTerm(
                METUtil.MuonTotal,
                event.MET_Muon_Total_Staco_etx,
                event.MET_Muon_Total_Staco_ety,
                event.MET_Muon_Total_Staco_sumet)

        # Note that RefMuon is not rebuilt from muons
        # -- it is a calorimeter term.
        self.met_utility.setMETTerm(
            METUtil.RefMuon,
            event.MET_RefMuon_Staco_etx,
            event.MET_RefMuon_Staco_ety,
            event.MET_RefMuon_Staco_sumet)

        """
        TAUS
        """
        if self.terms & Systematics.TAU_TERMS:
            self.met_utility.setTauParameters(
                event.tau_pt,
                event.tau_eta,
                event.tau_phi,
                event.tau_MET_wet,
                event.tau_MET_wpx,
                event.tau_MET_wpy,
                event.tau_MET_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefTau,
                event.MET_RefTau_etx,
                event.MET_RefTau_ety,
                event.MET_RefTau_sumet)

        self.met_utility.setMETTerm(
            METUtil.SoftTerms,
            event.MET_CellOut_Eflow_STVF_etx,
            event.MET_CellOut_Eflow_STVF_ety,
            event.MET_CellOut_Eflow_STVF_sumet)

        MET = self.get_met()

        if self.verbose:
            log.info("Run: {0} Event: {1}".format(
                event.RunNumber,
                event.EventNumber))
            log.info("Recalculated MET: %.3f (original: %.3f)" % (
                     MET.et(), event.MET_RefFinal_STVF_et))
            log.info("Recalculated MET phi: %.3f (original: %.3f)" % (
                     MET.phi(), event.MET_RefFinal_STVF_phi))
            if (abs(MET.et() - event.MET_RefFinal_STVF_et) /
                    event.MET_RefFinal_STVF_et) > 0.1:
                log.warning("Large MET difference!")

        # update the MET with the shifted value
        self.tree.MET_etx_original = event.MET_RefFinal_STVF_etx
        self.tree.MET_ety_original = event.MET_RefFinal_STVF_ety
        self.tree.MET_et_original = event.MET_RefFinal_STVF_et
        self.tree.MET_sumet_original = event.MET_RefFinal_STVF_sumet
        self.tree.MET_phi_original = event.MET_RefFinal_STVF_phi

        event.MET_RefFinal_STVF_etx = MET.etx()
        event.MET_RefFinal_STVF_ety = MET.ety()
        event.MET_RefFinal_STVF_et = MET.et()
        event.MET_RefFinal_STVF_sumet = MET.sumet()
        event.MET_RefFinal_STVF_phi = MET.phi()

        return True

    """
    Muon and photon systematics need to become their own classes like
    JES, TES above
    """

    def photon_systematics(self, event):
        """
        PHOTON SYSTEMATICS
        """
        # Now we get the same for photons
        self.pesUp.clear()
        self.pesDown.clear()
        self.perUp.clear()
        self.perDown.clear()
        self.ph_smeared_pt.clear()

        for iph, ph in enumerate(event.photons):

            self.egammaTool.SetRandomSeed(int(1e5*abs(ph.phi)))

            # Smear to match the data resolution, or by systematic variations
            smear = self.egammaTool.getSmearingCorrectionMeV(ph.cl_eta, ph.E, 0, True)
            smearUp = self.egammaTool.getSmearingCorrectionMeV(ph.cl_eta, ph.E, 2, True)
            smearDown = self.egammaTool.getSmearingCorrectionMeV(ph.cl_eta, ph.E, 1, True)

            self.ph_smeared_pt.push_back(smear * ph.pt)
            self.perUp.push_back((smearUp - smear) / smear)
            self.perDown.push_back((smearDown - smear) / smear)

            # Correct the measured energies in data, and scale by systematic variations
            # Conversions are treated differently.
            correction = 1.
            photontype = "CONVERTED_PHOTON" if ph.isConv else "UNCONVERTED_PHOTON"
            if self.datatype == datasets.DATA:
                correction = self.egammaTool.applyEnergyCorrectionMeV(
                        ph.cl_eta,
                        ph.cl_phi,
                        ph.E,
                        ph.cl_pt,
                        0,photontype) / ph.E

            self.ph_smeared_pt[iph] *= correction
            energyUp = self.egammaTool.applyEnergyCorrectionMeV(
                    ph.cl_eta,
                    ph.cl_phi,
                    ph.E,
                    ph.cl_pt,
                    2,photontype) / (correction * ph.E) - 1
            energyDown = self.egammaTool.applyEnergyCorrectionMeV(
                    ph.cl_eta,
                    ph.cl_phi,
                    ph.E,
                    ph.cl_pt,
                    1,photontype) / (correction * ph.E) - 1

            self.pesUp.push_back(energyUp)
            self.pesDown.push_back(energyDown)


    def muon_systematics(self, event):
        """
        MUON SYSTEMATICS
        """
        # And now the same for muons. We need resolution shifts for ID and MS,
        # and different treatment for the MS four-vector (for standalone muons).
        self.mu_smeared_pt.clear()
        self.mu_smeared_ms_pt.clear()

        self.cb_meridUp.clear()
        self.cb_meridDown.clear()
        self.cb_mermsUp.clear()
        self.cb_mermsDown.clear()
        self.mermsUp.clear()
        self.mermsDown.clear()

        self.mesUp.clear()
        self.mesDown.clear()

        for mu in event.muons:

            ptcb = mu.pt
            ptid = abs(sin(mu.id_theta_exPV)/mu.id_qoverp_exPV) if mu.id_qoverp_exPV != 0. else 0.
            ptms = abs(sin(mu.ms_theta)/mu.ms_qoverp) if mu.ms_qoverp != 0. else 0.
            self.muonTool.SetSeed(int(1.e+5*fabs(mu.phi)))
            etaMu = mu.eta
            charge = mu.charge
            self.muonTool.Event(ptms, ptid, ptcb, etaMu, charge)

            smearedCombinedPt = self.muonTool.pTCB()
            if not mu.isCombinedMuon:
                smearedCombinedPt = self.muonTool.pTMS() + self.muonTool.pTID()

            smearedMSPt = self.muonTool.pTMS()

            self.mu_smeared_ms_pt.push_back(smearedMSPt)
            self.mu_smeared_pt.push_back(smearedCombinedPt)

            smearedpTMS = 0.1
            smearedpTID = 0.1
            smearedpTCB = 0.1

            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSLOW")
            smearedpTMS = ptMS_smeared/smearedMSPt - 1.0
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            self.mermsDown.push_back(smearedpTMS)
            self.cb_mermsDown.push_back(smearedpTCB)
            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSUP")
            smearedpTMS = ptMS_smeared/smearedMSPt - 1.0
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            self.mermsUp.push_back(smearedpTMS)
            self.cb_mermsUp.push_back(smearedpTCB)
            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDUP")
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            self.cb_meridUp.push_back(smearedpTCB)
            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDLOW")
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            self.cb_meridDown.push_back(smearedpTCB)

            detRegion = self.muonTool.DetRegion()
            if detRegion == -1:
                detRegion = 3
            scalesyst = self.muonTool.getScaleSyst_CB().at(detRegion)
            self.mesUp.push_back(scalesyst)
            self.mesDown.push_back(-scalesyst)

        # Muons are complicated, and MET makes use of track, spectro, and combined quantites,
        # so we need all of their resolutions.
        # comboms reflects that it is the combined muon res affected by shifting ms res up and down.
        # comboid is for shifting the id res up and down
        self.met_utility.setObjectResolutionShift(
                METUtil.MuonsComboMS,
                self.cb_mermsUp,
                self.cb_mermsDown)

        self.met_utility.setObjectResolutionShift(
                METUtil.MuonsComboID,
                self.cb_meridUp,
                self.cb_meridDown)

        self.met_utility.setObjectResolutionShift(
                METUtil.SpectroMuons,
                self.mermsUp,
                self.mermsDown)

        # For now the mes only affects combined
        self.met_utility.setObjectEnergyUncertainties(
                METUtil.Muons,
                self.mesUp,
                self.mesDown)

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
from . import utils
from . import datasets
from .units import GeV
from .rand import get_random

# rootpy imports
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
from externaltools import JVFUncertaintyTool as JVFUncertaintyTool2012

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
        self.tree = sys_util.tree
        self.verbose = sys_util.verbose
        self.year = sys_util.year


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
        return [ self.jes_tool.getRelPileupPtTerm(jet.pt, jet.eta, self.tree.nvtxjets, event.averageIntPerXing) ]

    def _PUNPV(self, jet, event):
        return [ self.jes_tool.getRelNPVOffsetTerm(jet.pt, jet.eta, self.tree.nvtxjets) ]

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
        # JVF is only used in a certain range, so only correct for those
        if jet.pt < 50e3 and abs(jet.constscale_eta) < 2.4:
            truejets = VectorTLorentzVector()
            truejets_cache = []
            for truejet in event.truejets:
                if truejet.pt > 10e3:
                    t = TLorentzVector()
                    t.SetPtEtaPhiM(truejet.pt, truejet.eta,
                                   truejet.phi, truejet.m)
                    truejets.push_back(t)
                    truejets_cache.append(t)
            isPU = self.jvf_tool.isPileUpJet(jet.fourvect, truejets)
            jvf_cut_sys = self.jvf_tool.getJVFcut(
                self.JVFcutNominal, isPU,
                jet.pt, jet.constscale_eta, self.is_up)
            jvf_cut_diff = jvf_cut_sys - self.JVFcutNominal
            jet.jvtxf -= jvf_cut_diff


class JER(JetSystematic):

    def __init__(self, is_up, **kwargs):
        super(JER, self).__init__(is_up, **kwargs)
        # Tag assumed: JetResolution-01-00-00
        if self.year == 2011:
            log.info("Using 2011 JER config")
            self.jer_tool = JERProvider(
                "AntiKt4LCTopoJES", "Truth",
                JetResolution.get_resource('JERProviderPlots_2011.root'))
            self.jer_tool.is7TeV(True)
        else:
            log.info("Using 2012 JER config")
            self.jer_tool = JERProvider(
                "AntiKt4LCTopoJES", "Truth",
                JetResolution.get_resource('JERProviderPlots_2012.root'))
        self.jer_tool.init()
        self.jetrandom = get_random()

    @JetSystematic.set
    def run(self, jet, event):
        """
        Note: The JERDown shift is essentially meaningless.
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/JetEnergyResolutionProvider2012
        """
        shift = 1.
        # Allowable range is > 10 GeV, but anything below 20 enters SoftJets
        if jet.pt > 20 * GeV and jet.pt < 10000 * GeV:
            smear = self.jer_tool.getSmearingFactorMC(jet.pt, jet.eta)
            shift = self.jetrandom.Gaus(1., smear)
        jet.pt *= shift


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
                 np='TOTAL',
                 matched_state=None,
                 **kwargs):
        super(TES, self).__init__(is_up, **kwargs)
        self.np = getattr(TauCorrUncert.TESUncertainty, np)
        self.matched_state = matched_state
        if self.year == 2011:
            input = TCU.get_resource('TES/mc11.root')
        elif self.year == 2012:
            # TODO use medium and tight?
            input = TCU.get_resource('TES/mc12_p1344_medium.root')
        else:
            raise ValueError('No TES defined for year %d' % self.year)
        self.tes_tool = TauCorrUncert.TESUncertainty(input)

    @TauSystematic.set
    def run(self, tau, event):
        # need nominal pT for trigger efficiency correction
        tau._pt_nominal = tau.pt
        if self.matched_state is True and not tau.matched:
            return
        elif self.matched_state is False and tau.matched:
            return
        shift = self.tes_tool.GetTESUncertainty(
            tau.pt, tau.eta, tau.numTrack, self.np)
        if self.is_up:
            tau.pt *= 1. + shift
        else:
            tau.pt *= 1. - shift


class Systematics(EventFilter):

    # default
    NONE = METUtil.None

    # taus
    # TES 2011
    TES_TRUE_FINAL_UP = -30000
    TES_TRUE_FINAL_DOWN = -30001
    TES_FAKE_FINAL_UP = -30002
    TES_FAKE_FINAL_DOWN = -30003
    # TES 2012
    TES_TRUE_INSITUINTERPOL_UP = -30004
    TES_TRUE_INSITUINTERPOL_DOWN = -30005
    TES_TRUE_SINGLEPARTICLEINTERPOL_UP = -30006
    TES_TRUE_SINGLEPARTICLEINTERPOL_DOWN = -30007
    TES_TRUE_MODELING_UP = -30008
    TES_TRUE_MODELING_DOWN = -30009
    TES_FAKE_TOTAL_UP = -30010
    TES_FAKE_TOTAL_DOWN = -30011

    TES_TERMS = set([
        TES_TRUE_FINAL_UP,
        TES_TRUE_FINAL_DOWN,
        TES_FAKE_FINAL_UP,
        TES_FAKE_FINAL_DOWN,
        TES_TRUE_INSITUINTERPOL_UP,
        TES_TRUE_INSITUINTERPOL_DOWN,
        TES_TRUE_SINGLEPARTICLEINTERPOL_UP,
        TES_TRUE_SINGLEPARTICLEINTERPOL_DOWN,
        TES_TRUE_MODELING_UP,
        TES_TRUE_MODELING_DOWN,
        TES_FAKE_TOTAL_UP,
        TES_FAKE_TOTAL_DOWN,
    ])

    TER_UP = METUtil.TERUp
    TER_DOWN = METUtil.TERDown

    TAU_TERMS = set([TER_UP, TER_DOWN]) | TES_TERMS

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

    def __init__(self,
            datatype,
            year,
            tree,
            terms=None,
            verbose=False,
            **kwargs):

        super(Systematics, self).__init__(**kwargs)

        self.systematics = []
        self.terms = set([])
        self.datatype = datatype
        self.year = year
        self.tree = tree
        self.verbose = verbose

        if terms is not None:
            # remove possible duplicates
            terms = set(terms)
            self.terms = terms
            for term in terms:
                systematic = None

                # JES terms
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

                # JVF
                elif term == Systematics.JVF_UP:
                    systematic = JVF(True, sys_util=self)
                elif term == Systematics.JVF_DOWN:
                    systematic = JVF(False, sys_util=self)

                # JER (one-sided)
                elif term == Systematics.JER_UP:
                    systematic = JER(True, sys_util=self)

                # TES 2011 true up and down
                elif term == Systematics.TES_TRUE_FINAL_UP:
                    systematic = TES(True, sys_util=self,
                        np='FINAL',
                        matched_state=True)
                elif term == Systematics.TES_TRUE_FINAL_DOWN:
                    systematic = TES(False, sys_util=self,
                        np='FINAL',
                        matched_state=True)

                # TES 2012 true up and down
                elif term == Systematics.TES_TRUE_INSITUINTERPOL_UP:
                    systematic = TES(True, sys_util=self,
                        np='INSITUINTERPOL',
                        matched_state=True)
                elif term == Systematics.TES_TRUE_INSITUINTERPOL_DOWN:
                    systematic = TES(False, sys_util=self,
                        np='INSITUINTERPOL',
                        matched_state=True)
                elif term == Systematics.TES_TRUE_SINGLEPARTICLEINTERPOL_UP:
                    systematic = TES(True, sys_util=self,
                        np='SINGLEPARTICLEINTERPOL',
                        matched_state=True)
                elif term == Systematics.TES_TRUE_SINGLEPARTICLEINTERPOL_DOWN:
                    systematic = TES(False, sys_util=self,
                        np='SINGLEPARTICLEINTERPOL',
                        matched_state=True)
                elif term == Systematics.TES_TRUE_MODELING_UP:
                    systematic = TES(True, sys_util=self,
                        np='MODELING',
                        matched_state=True)
                elif term == Systematics.TES_TRUE_MODELING_DOWN:
                    systematic = TES(False, sys_util=self,
                        np='MODELING',
                        matched_state=True)

                # TES 2011 fake up and down
                elif term == Systematics.TES_FAKE_FINAL_UP:
                    systematic = TES(True, sys_util=self,
                        np='FINAL',
                        matched_state=False)
                elif term == Systematics.TES_FAKE_FINAL_DOWN:
                    systematic = TES(False, sys_util=self,
                        np='FINAL',
                        matched_state=False)

                # TES 2012 fake up and down
                elif term == Systematics.TES_FAKE_TOTAL_UP:
                    systematic = TES(True, sys_util=self,
                        np='TOTAL',
                        matched_state=False)
                elif term == Systematics.TES_FAKE_TOTAL_DOWN:
                    systematic = TES(False, sys_util=self,
                        np='TOTAL',
                        matched_state=False)

                # MET terms handled in met.py
                elif term not in Systematics.MET_TERMS:
                    raise ValueError("systematic not supported")

                if systematic is not None:
                    self.systematics.append(systematic)

    def passes(self, event):
        # run the systematics
        for systematic in self.systematics:
            systematic.run(event)
        return True

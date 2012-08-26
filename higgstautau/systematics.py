"""
Adapted from the Example.C in MissingETUtility/macros
"""
# stdlib Python imports
import math
from math import sin, sqrt, pow

# local imports
from . import tauid

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

# MissingETUtility
from ROOT import METUtility
from ROOT import METUtil
from ROOT import METObject
# JetUncertainties
from ROOT import MultijetJESUncertaintyProvider
# JetResolution
from ROOT import JERProvider
# egammaAnalysisUtils
from ROOT import eg2011
# MuonMomentumCorrections
from ROOT import MuonSmear
# MissingETUtility
from ROOT import TESUncertaintyProvider


class ObjectSystematic(object):

    def __init__(self, sys_util, verbose=False):

        self.sys_util = sys_util
        self.verbose = verbose


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
                    print jet.pt
                print "-" * 20

            for jet in event.jets:
                f(self, jet, event)

            if self.verbose:
                print "JETS AFTER:"
                for jet in event.jets:
                    print jet.pt

        return wrapper


class JES(JetSystematic):

    def __init__(self, is_up, **kwargs):

        # *** Set up the uncertainty tools ***
        # Tag assumed: JetUncertainties-00-05-09-02
        self.jes_tool = MultijetJESUncertaintyProvider(
            JetUncertainties.get_resource("MultijetJES_Preliminary.config"),
            JetUncertainties.get_resource("InsituJES2011_AllNuisanceParameters.config"),
            "AntiKt4LCTopoJets","MC11c")
        super(JES, self).__init__(is_up, **kwargs)

    @JetSystematic.set
    def run(self, jet, event):

        # Safest to assume nothing about the uncertainties on soft jets.
        # These will go into SoftJets anyhow, and so the JES systematics
        # aren't used.
        shift = 0.
        if jet.pt > 20e3 and jet.pt < 7000e3 and abs(jet.eta) < 4.5:

            # delta R cut needed to apply close-by jets uncertainty
            drmin=9999
            for otherjet in event.jets:
                if otherjet.emscale_pt > 7000:
                    if jet.index != otherjet.index:
                        dr = utils.dR(jet.eta, jet.phi,
                                      otherjet.eta,
                                      otherjet.phi)
                        if dr < drmin:
                            drmin = dr

            # The bool is the "isPos" argument
            if self.is_up:
                shift = self.jes_tool.getRelUncert(jet.pt,
                             jet.eta, drmin,
                             True, self.sys_util.nvtxjets, event.averageIntPerXing)
            else:
                shift = -1 * self.jes_tool.getRelUncert(jet.pt,
                              jet.eta, drmin,
                              False, self.sys_util.nvtxjets, event.averageIntPerXing)
        jet.pt *= 1. + shift


class JER(JetSystematic):

    def __init__(self, is_up, **kwargs):

        # Tag assumed: JetResolution-01-00-00
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

    def __init__(self, is_up, **kwargs):

        # No tag yet, testing code
        self.tes_tool = TESUncertaintyProvider()
        super(TES, self).__init__(is_up, **kwargs)

    @TauSystematic.set
    def run(self, tau, event):

        pt = tau.pt
        eta = tau.eta
        nProng = tau.nProng
        shift = self.tes_tool.GetTESUncertainty(pt / 1e3, eta, nProng)
        if shift < 0:
            shift = 0
        if not self.is_up:
            shift *= -1
        # need nominal pT for trigger efficiency correction
        tau._pt_nominal = tau.pt
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
    ELECTRON_TERMS = {EES_UP, EES_DOWN, EER_UP, EER_DOWN}

    # photons
    PES_UP = METUtil.PESUp
    PES_DOWN = METUtil.PESDown
    PER_UP = METUtil.PERUp
    PER_DOWN = METUtil.PERDown
    PHOTON_TERMS = {PES_UP, PES_DOWN, PER_UP, PER_DOWN}

    # taus
    TES_UP = METUtil.TESUp
    TES_DOWN = METUtil.TESDown
    TER_UP = METUtil.TERUp
    TER_DOWN = METUtil.TERDown
    TES_TERMS = {TES_UP, TES_DOWN}
    TAU_TERMS = {TES_UP, TES_DOWN, TER_UP, TER_DOWN}
    TAUBDT_UP = -100
    TAUBDT_DOWN = -101

    # jets
    JES_UP = METUtil.JESUp
    JES_DOWN = METUtil.JESDown
    JER_UP = METUtil.JERUp
    JER_DOWN = METUtil.JERDown # NOT USED!
    JET_TERMS = {JES_UP, JES_DOWN, JER_UP, JER_DOWN}

    # muons
    MERID_UP = METUtil.MERIDUp
    MERID_DOWN = METUtil.MERIDDown
    MERMS_UP = METUtil.MERMSUp
    MERMS_DOWN = METUtil.MERMSDown
    MES_UP = METUtil.MESUp
    MES_DOWN = METUtil.MESDown
    MUON_TERMS = {MERID_UP, MERID_DOWN, MERMS_UP, MERMS_DOWN, MES_UP, MES_DOWN}

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

    ALLCLUSTERS_UP = METUtil.AllClustersUp
    ALLCLUSTERS_DOWN = METUtil.AllClustersDown

    # soft terms
    RESOSOFTTERMS_PTHARD_UP = METUtil.ResoSoftTermsUp_ptHard
    RESOSOFTTERMS_PTHARD_DOWN = METUtil.ResoSoftTermsDown_ptHard

    RESOSOFTTERMS_PTHARD_UPDOWN = METUtil.ResoSoftTermsUpDown_ptHard
    RESOSOFTTERMS_PTHARD_DOWNUP = METUtil.ResoSoftTermsDownUp_ptHard

    SCALESOFTTERMS_PTHARD_UP = METUtil.ScaleSoftTermsUp_ptHard
    SCALESOFTTERMS_PTHARD_DOWN = METUtil.ScaleSoftTermsDown_ptHard

    RESOSOFTTERMS_UP = METUtil.ResoSoftTermsUp
    RESOSOFTTERMS_DOWN = METUtil.ResoSoftTermsDown

    SCALESOFTTERMS_UP = METUtil.ScaleSoftTermsUp
    SCALESOFTTERMS_DOWN = METUtil.ScaleSoftTermsDown

    # pileup
    PILEUP_UP = METUtil.PileupUp
    PILEUP_DOWN = METUtil.PileupDown

    def __init__(self,
            datatype,
            year,
            terms=None,
            verbose=False,
            very_verbose=False,
            **kwargs):

        super(Systematics, self).__init__(**kwargs)
        self.systematics = []
        self.terms = set([])

        if terms is not None:
            # remove possible duplicates
            terms = set(terms)
            self.terms = terms
            for term in terms:
                if term == Systematics.JES_UP:
                    systematic = JES(True, sys_util=self, verbose=verbose)
                elif term == Systematics.JES_DOWN:
                    systematic = JES(False, sys_util=self, verbose=verbose)
                elif term == Systematics.JER_UP:
                    systematic = JER(True, sys_util=self, verbose=verbose)
                elif term == Systematics.TES_UP:
                    systematic = TES(True, sys_util=self, verbose=verbose)
                elif term == Systematics.TES_DOWN:
                    systematic = TES(False, sys_util=self, verbose=verbose)
                elif term == Systematics.EES_UP:
                    systematic = EES(True, datatype=datatype, sys_util=self,
                            verbose=verbose)
                elif term == Systematics.EES_DOWN:
                    systematic = EES(False, datatype=datatype, sys_util=self,
                            verbose=verbose)
                elif term == Systematics.EER_UP:
                    systematic = EER(True, sys_util=self, verbose=verbose)
                elif term == Systematics.EER_DOWN:
                    systematic = EER(False, sys_util=self, verbose=verbose)
                elif term == Systematics.TAUBDT_UP:
                    systematic = TauBDT(True, sys_util=self, verbose=verbose)
                elif term == Systematics.TAUBDT_DOWN:
                    systematic = TauBDT(False, sys_util=self, verbose=verbose)
                else:
                    raise ValueError("systematic not supported")
                self.systematics.append(systematic)

        self.datatype = datatype
        self.year = year
        self.verbose = verbose
        self.very_verbose = very_verbose

        # Initialise your METUtility object
        self.met_utility = METUtility()

        # *** Demonstration of the configuration here  ***
        # *** All values that are set are the defaults ***

        # Turn on (off) the relevant MET terms
        # These are the terms required for MET_RefFinal(_BDTMedium)
        self.met_utility.defineMissingET(
                True,  # RefEle
                True,  # RefGamma
                True,  # RefTau
                True,  # RefJet
                True,  # SoftJets
                True,  # RefMuon
                True,  # MuonTotal
                False, # CellOut
                True   # CellOut_Eflow
            )

        # The threshold below which jets enter the SoftJets term (JES is not applied)
        self.met_utility.setSoftJetCut(20e3)

        # Whether to use MUID muons (otherwise STACO).
        self.met_utility.setIsMuid(False)

        # Whether METUtility should scream at you over every little thing
        self.met_utility.setVerbosity(self.very_verbose)

        # Tag assumed: MuonMomentumCorrections-00-05-03
        self.muonTool = MuonSmear.SmearingClass(
                "Data11","staco","pT","Rel17",
                MuonMomentumCorrections.RESOURCE_PATH)
        self.muonTool.UseScale(1)
        self.muonTool.UseImprovedCombine()

    def passes(self, event):
        """
        Demonstrates how to set up the METUtility with object momenta
        such that MET_RefFinal can be rebuilt matching the values in D3PD.

        *** *** *** *** *** DISCLAIMER *** *** *** *** ***

        These examples of uncertainty-setting are meant to
        demonstrate how the uncertainties should be passed
        to METUtility.  Recommendations on just which ones
        you are meant to be using come from the CP groups.
        """
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

        """
        JETS

        Always use setJetParameters since they may be recalibrated upstream
        """
        self.met_utility.setJetParameters(
            event.jet_pt,
            event.jet_eta,
            event.jet_phi,
            event.jet_E,
            event.jet_AntiKt4LCTopo_MET_BDTMedium_wet,
            event.jet_AntiKt4LCTopo_MET_BDTMedium_wpx,
            event.jet_AntiKt4LCTopo_MET_BDTMedium_wpy,
            event.jet_AntiKt4LCTopo_MET_BDTMedium_statusWord)

        self.met_utility.setOriJetParameters(event.jet_pt)

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

        self.met_utility.setMETTerm(
            METUtil.SoftJets,
            event.MET_SoftJets_BDTMedium_etx,
            event.MET_SoftJets_BDTMedium_ety,
            event.MET_SoftJets_BDTMedium_sumet)

        self.met_utility.setMETTerm(
            METUtil.CellOutEflow,
            event.MET_CellOut_BDTMedium_etx,
            event.MET_CellOut_BDTMedium_ety,
            event.MET_CellOut_BDTMedium_sumet)

        MET = self.met_utility.getMissingET(METUtil.RefFinal)

        if self.verbose:
            print "Recalculated MET: %.3f (original: %.3f)" % (
                    MET.et(), event.MET_RefFinal_BDTMedium_et)

        # update the MET with the shifted value
        event.MET_RefFinal_BDTMedium_etx = MET.etx()
        event.MET_RefFinal_BDTMedium_ety = MET.ety()
        event.MET_RefFinal_BDTMedium_et = MET.et()
        event.MET_RefFinal_BDTMedium_sumet = MET.sumet()
        event.MET_RefFinal_BDTMedium_phi = MET.phi()

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


    def consistency(self, event):

        # Demonstrates how to set up the METUtility with object momenta such that
        # MET_RefFinal can be rebuilt matching the values in D3PD.

        # Start with a clean METUtility
        self.testUtil.reset()

        # For these, just the kinematics need to be set
        self.testUtil.setElectronParameters(
                event.el_pt,
                event.el_eta,
                event.el_phi,
                event.el_MET_wet,
                event.el_MET_wpx,
                event.el_MET_wpy,
                event.el_MET_statusWord)

        self.testUtil.setPhotonParameters(
                event.ph_pt,
                event.ph_eta,
                event.ph_phi,
                event.ph_MET_wet,
                event.ph_MET_wpx,
                event.ph_MET_wpy,
                event.ph_MET_statusWord)

        #  Tau rebuilding is unsafe. The tau finding is frequently rerun in D3PD.
        #  self.testUtil.setTauParameters(tau_pt, tau_eta, tau_phi,
        #			     tau_MET_wet, tau_MET_wpx, tau_MET_wpy,
        #			     tau_MET_statusWord)

        #  Cluster rebuilding is just for specialist studies
        #  self.testUtil.setClusterParameters(cl_pt, cl_eta, cl_phi,
        #				 cl_MET_wet, cl_MET_wpx, cl_MET_wpy,
        #				 cl_MET_statusWord)

        #  Track rebuilding is just for specialist studies
        #  self.testUtil.setTrackParameters(trk_pt, trk_eta, trk_phi,
        #			       trk_MET_wet, trk_MET_wpx, trk_MET_wpy,
        #			       trk_MET_statusWord)

        # The SoftJets term is now to be taken from D3PD, and no "hard jets" are allowed
        # to enter it. Recalibration or smearing could cause jets that were above the
        # 20 GeV threshold to drop below it, so we supply the original pT's to prevent
        # them from being moved out of RefJet.
        self.testUtil.setJetParameters(
                event.jet_pt,
                event.jet_eta,
                event.jet_phi,
                event.jet_E,
                event.jet_MET_wet,
                event.jet_MET_wpx,
                event.jet_MET_wpy,
                event.jet_MET_statusWord)
        self.testUtil.setOriJetParameters(event.jet_pt)

        # Muons may be ID, combined, or standalone. For the latter especially,
        # we need to set the MS four-momenta because they are treated differently.
        self.testUtil.setMuonParameters(
                event.mu_pt,
                event.mu_eta,
                event.mu_phi,
                event.mu_MET_wet,
                event.mu_MET_wpx,
                event.mu_MET_wpy,
                event.mu_MET_statusWord)
        self.testUtil.setExtraMuonParameters(
                event.mu_ms_qoverp,
                event.mu_ms_theta,
                event.mu_ms_phi,
                event.mu_charge)
        # An alternative version of this method is demonstrated below, and takes pT, eta, phi instead.
        # This is more convenient when one needs to smear the pT, for example.

        # When the terms are not rebuilt from the objects, due to incomplete info,
        # then one needs to set the term directly from the stored value in the D3PD.
        # This might also be done if you aren't interested in the possible variations
        # of that term. E.g. if you only care about photons, no need to rebuild muons.

        # self.testUtil.setMETTerm(METUtil.RefJet, MET_RefJet_etx, MET_RefJet_ety, MET_RefJet_sumet)
        self.testUtil.setMETTerm(
                METUtil.SoftJets,
                event.MET_SoftJets_etx,
                event.MET_SoftJets_ety,
                event.MET_SoftJets_sumet)
        # self.testUtil.setMETTerm(METUtil.RefEle, MET_RefEle_etx, MET_RefEle_ety, MET_RefEle_sumet)
        # self.testUtil.setMETTerm(METUtil.RefGamma, MET_RefGamma_etx, MET_RefGamma_ety, MET_RefGamma_sumet)

        # *** Note the difference in naming -- there is normally no MET_MuonTotal term.
        #     It's usually Muon_Total, Muon_Total_Muid, something like that.
        #     MET_RefFinal in particular uses MET_MuonBoy.
        # self.testUtil.setMETTerm(METUtil.MuonTotal, MET_MuonBoy_etx, MET_MuonBoy_ety, MET_MuonBoy_sumet)

        # *** Note that RefMuon is not rebuilt from muons -- it is a calorimeter term.
        self.testUtil.setMETTerm(
                METUtil.RefMuon,
                event.MET_RefMuon_etx,
                event.MET_RefMuon_ety,
                event.MET_RefMuon_sumet)
        self.testUtil.setMETTerm(
                METUtil.RefTau,
                event.MET_RefTau_etx,
                event.MET_RefTau_ety,
                event.MET_RefTau_sumet)
        # self.testUtil.setMETTerm(METUtil.CellOut, MET_CellOut_etx, MET_CellOut_ety, MET_CellOut_sumet)
        self.testUtil.setMETTerm(
                METUtil.CellOutEflow,
                event.MET_CellOut_Eflow_etx,
                event.MET_CellOut_Eflow_ety,
                event.MET_CellOut_Eflow_sumet)

        # This is the simple check, where you compare the terms manually against what's in the D3PD.
        # Note: every call to getMissingET recomputes the terms, so if you need to get more than one
        # value, e.g. etx, ety, et, sumET, it's more efficient to get the METObject.
        # Usually, comparing etx and/or ety is more informative, because et could be right if
        # etx and ety were flipped, for example. They also add linearly, which et doesn't.

        RefEle_util = self.testUtil.getMissingET(METUtil.RefEle)
        RefGamma_util = self.testUtil.getMissingET(METUtil.RefGamma)
        RefTau_util = self.testUtil.getMissingET(METUtil.RefTau)
        RefMuon_util = self.testUtil.getMissingET(METUtil.RefMuon)
        RefJet_util = self.testUtil.getMissingET(METUtil.RefJet)
        SoftJets_util = self.testUtil.getMissingET(METUtil.SoftJets)
        #refCellOut_util = self.testUtil.getMissingET(METUtil.CellOut)
        CellOutEflow_util = self.testUtil.getMissingET(METUtil.CellOutEflow)
        MuonTotal_util = self.testUtil.getMissingET(METUtil.MuonTotal)
        RefFinal_util = self.testUtil.getMissingET(METUtil.RefFinal)

        if self.verbose:
            print >> self.stream, "** Manual consistency check **" << endl
            print >> self.stream, "Term:    Original   vs    Tool output"
            print >> self.stream, "RefEle etx: "    << MET_RefEle_etx   << " vs " << RefEle_util.etx()
            print >> self.stream, "RefGamma etx: "  << MET_RefGamma_etx << " vs " << RefGamma_util.etx()
            print >> self.stream, "RefTau etx: "    << MET_RefTau_etx   << " vs " << RefTau_util.etx()
            print >> self.stream, "RefMuon etx: "   << MET_RefMuon_etx  << " vs " << RefMuon_util.etx()
            print >> self.stream, "RefJet etx: "    << MET_RefJet_etx   << " vs " << RefJet_util.etx()
            print >> self.stream, "SoftJets etx: "   << MET_SoftJets_etx << " vs " << SoftJets_util.etx()
            print >> self.stream, "MuonBoy etx: "   << MET_MuonBoy_etx  << " vs " << MuonTotal_util.etx()
            print >> self.stream, "CellOut_Eflow etx: " << MET_CellOut_Eflow_etx << " vs " << CellOutEflow_util.etx()
            print >> self.stream, "RefFinal etx: "  << MET_RefFinal_etx << " vs " << RefFinal_util.etx()  << endl

        # If you don't want to test manually, there's a built-in consistency check.
        # To test just one term, fill a METObject with the appropriate values,
        # then feed it to the checkConsistency() method.
        # The difference can be retrieved via a reference argument.

        refFinal_test = METObject(MET_RefFinal_etx,
                  MET_RefFinal_ety,
                  MET_RefFinal_sumet)
        check_refFinal = self.testUtil.checkConsistency(refFinal_test,METUtil.RefFinal)
        if check_refFinal:
            print >> self.stream, "RefFinal checks out!"
        else:
            print >> self.stream, "RefFinal doesn't check out!"

        # By filling a vector of terms, you can test all of them in one function call.
        # The sum (and sumET) will be tested as well. (Can't get the difference this way).
        refEle_test = METObject(MET_RefEle_etx,
                MET_RefEle_ety,
                MET_RefEle_sumet)
        refGamma_test = METObject(MET_RefGamma_etx,
                MET_RefGamma_ety,
                MET_RefGamma_sumet)
        refJet_test = METObject(MET_RefJet_etx,
                MET_RefJet_ety,
                MET_RefJet_sumet)
        muonBoy_test = METObject(MET_MuonBoy_etx,
                 MET_MuonBoy_ety,
                 MET_MuonBoy_sumet)

        testvector = ROOT.vector('pair<int,METObject>')()
        testvector.push_back(ROOT.pair('int,METObject')(METUtil.RefEle,refEle_test))
        testvector.push_back(ROOT.pair('int,METObject')(METUtil.RefGamma,refGamma_test))
        testvector.push_back(ROOT.pair('int,METObject')(METUtil.RefJet,refJet_test))
        testvector.push_back(ROOT.pair('int,METObject')(METUtil.MuonTotal,muonBoy_test))

        check = self.testUtil.checkConsistency(testvector)
        if check:
            print >> self.stream, "MET checks out!"
        else:
            print >> self.stream, "MET doesn't check out!"

        # In addition to the etx, ety, sumet retrieval, you can also get the MET significance.
        # By default, METObject.sig() returns etx() / ( 0.5*sqrt(sumet()) )
        #
        # There is also the possibility of returning a more sophisticated estimator for
        # the significance, activated by calling METUtility.doSignificance() in setup.
        # This is still under development, and requires all object resolutions to be set.


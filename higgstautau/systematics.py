"""
Adapted from the Example.C in MissingETUtility/macros
"""
# stdlib Python imports
import os
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
# load TESUncertaintyProvider before MissingETUtility!!!
from externaltools import TESUncertaintyProvider as TESP
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
# TESUncertaintyProvider
from ROOT import TESUncertaintyProvider


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
        super(JES, self).__init__(is_up, **kwargs)

        self.jes_tool = None

        if self.year == 2011:
            self.jes_tool = MultijetJESUncertaintyProvider(
                "JES_2011/Final/MultijetJES_2011.config",
                "JES_2011/Final/InsituJES2011_AllNuisanceParameters.config",
                "AntiKt4LCTopoJets",
                "MC11c",
                JetUncertainties.RESOURCE_PATH)
        elif self.year == 2012:
            self.jes_tool = MultijetJESUncertaintyProvider(
                "JES_2012/MultijetJES_2012.config",
                "JES_2012/JESUncertainty2012_Sept2012.config",
                "AntiKt4LCTopoJets",
                "MC12a",
                JetUncertainties.RESOURCE_PATH)
        else:
            raise ValueError('No JES defined for year %d' % year)

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

            # TODO: shift is symmetric (is_up argument is not needed)
            if self.is_up:
                shift = self.jes_tool.getRelUncert(
                        jet.pt,
                        jet.eta,
                        drmin,
                        True, # is up
                        self.sys_util.nvtxjets,
                        event.averageIntPerXing,
                        False) # is b jet
            else:
                shift = -1 * self.jes_tool.getRelUncert(
                        jet.pt,
                        jet.eta,
                        drmin,
                        False, # is up
                        self.sys_util.nvtxjets,
                        event.averageIntPerXing,
                        False) # is b jet

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

    def __init__(self, is_up, **kwargs):

        super(TES, self).__init__(is_up, **kwargs)

        if self.year == 2011:
            self.tes_tool = TESUncertaintyProvider(
                    os.path.normpath(TESP.RESOURCE_PATH), '', 'mc11')
        elif self.year == 2012:
            self.tes_tool = TESUncertaintyProvider(
                    os.path.normpath(TESP.RESOURCE_PATH), '', 'mc12')
        else:
            raise ValueError('No TES defined for year %d' % self.year)

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
    TER_UP = METUtil.TERUp
    TER_DOWN = METUtil.TERDown
    TES_TERMS = set([TES_UP, TES_DOWN])
    TAU_TERMS = set([TES_UP, TES_DOWN, TER_UP, TER_DOWN])
    TAUBDT_UP = -100
    TAUBDT_DOWN = -101

    # jets
    JES_UP = METUtil.JESUp
    JES_DOWN = METUtil.JESDown
    JER_UP = METUtil.JERUp
    JER_DOWN = METUtil.JERDown # NOT USED!
    JES_TERMS = set([JES_UP, JES_DOWN, JER_UP, JER_DOWN])

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

    # pileup (deprecated 1.1.0)
    #PILEUP_UP = METUtil.PileupUp
    #PILEUP_DOWN = METUtil.PileupDown

    def __init__(self,
            datatype,
            year,
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
        self.channel = channel
        self.verbose = verbose
        self.very_verbose = very_verbose

        if terms is not None:
            # remove possible duplicates
            terms = set(terms)
            self.terms = terms
            for term in terms:
                if term == Systematics.JES_UP:
                    systematic = JES(True, sys_util=self)
                elif term == Systematics.JES_DOWN:
                    systematic = JES(False, sys_util=self)
                elif term == Systematics.JER_UP:
                    systematic = JER(True, sys_util=self)
                elif term == Systematics.TES_UP:
                    systematic = TES(True, sys_util=self)
                elif term == Systematics.TES_DOWN:
                    systematic = TES(False, sys_util=self)
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
                else:
                    raise ValueError("systematic not supported")
                self.systematics.append(systematic)

        # Initialise your METUtility object
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
        # TODO: set the proper year here...
        #self.muonTool = MuonSmear.SmearingClass(
        #        "Data11","staco","pT","Rel17",
        #        MuonMomentumCorrections.RESOURCE_PATH)
        #self.muonTool.UseScale(1)
        #self.muonTool.UseImprovedCombine()

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
        """

        if self.channel == 'hh':
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

            self.met_utility.setOriJetParameters(event.jet_EtaOriginEM)

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

    def passes_12(self, event):
        """
        JETS
        Always use setJetParameters since they may be recalibrated upstream
        """
        self.met_utility.setJetParameters(
            event.jet_pt,
            event.jet_eta,
            event.jet_phi,
            event.jet_E,
            event.jet_AntiKt4LCTopo_MET_STVF_wet,
            event.jet_AntiKt4LCTopo_MET_STVF_wpx,
            event.jet_AntiKt4LCTopo_MET_STVF_wpy,
            event.jet_AntiKt4LCTopo_MET_STVF_statusWord)

        self.met_utility.setOriJetParameters(event.jet_pt)

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
                event.el_MET_STVF_wet,
                event.el_MET_STVF_wpx,
                event.el_MET_STVF_wpy,
                event.el_MET_STVF_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefEle,
                event.MET_RefEle_STVF_etx,
                event.MET_RefEle_STVF_ety,
                event.MET_RefEle_STVF_sumet)


        if self.terms & Systematics.PHOTON_TERMS:
            self.met_utility.setPhotonParameters(
                event.ph_pt,
                event.ph_eta,
                event.ph_phi,
                event.ph_MET_STVF_wet,
                event.ph_MET_STVF_wpx,
                event.ph_MET_STVF_wpy,
                event.ph_MET_STVF_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefGamma,
                event.MET_RefGamma_STVF_etx,
                event.MET_RefGamma_STVF_ety,
                event.MET_RefGamma_STVF_sumet)

        """
        MUONS
        """
        if self.terms & Systematics.MUON_TERMS:
            self.met_utility.setMuonParameters(
                event.mu_staco_pt, # or smeared pT
                event.mu_staco_eta,
                event.mu_staco_phi,
                event.mu_staco_MET_STVF_wet,
                event.mu_staco_MET_STVF_wpx,
                event.mu_staco_MET_STVF_wpy,
                event.mu_staco_MET_STVF_statusWord)

            # In this instance there is an overloaded version of
            # setExtraMuonParameters that accepts smeared pTs for spectro
            self.met_utility.setExtraMuonParameters(
                event.mu_staco_ms_pt, # or smeared pT
                event.mu_staco_ms_theta,
                event.mu_staco_ms_phi)
        else:
            self.met_utility.setMETTerm(
                METUtil.MuonTotal,
                event.MET_Muon_Total_Staco_STVF_etx,
                event.MET_Muon_Total_Staco_STVF_ety,
                event.MET_Muon_Total_Staco_STVF_sumet)

        # Note that RefMuon is not rebuilt from muons
        # -- it is a calorimeter term.
        self.met_utility.setMETTerm(
            METUtil.RefMuon,
            event.MET_RefMuon_Staco_STVF_etx,
            event.MET_RefMuon_Staco_STVF_ety,
            event.MET_RefMuon_Staco_STVF_sumet)

        """
        TAUS
        """
        if self.terms & Systematics.TAU_TERMS:
            self.met_utility.setTauParameters(
                event.tau_pt,
                event.tau_eta,
                event.tau_phi,
                event.tau_MET_STVF_wet,
                event.tau_MET_STVF_wpx,
                event.tau_MET_STVF_wpy,
                event.tau_MET_STVF_statusWord)
        else:
            self.met_utility.setMETTerm(
                METUtil.RefTau,
                event.MET_RefTau_STVF_etx,
                event.MET_RefTau_STVF_ety,
                event.MET_RefTau_STVF_sumet)

        #self.met_utility.setMETTerm(
        #    METUtil.SoftJets,
        #    event.MET_SoftJets_STVF_etx,
        #    event.MET_SoftJets_STVF_ety,
        #    event.MET_SoftJets_STVF_sumet) NOT NEEDED??

        self.met_utility.setMETTerm(
            METUtil.CellOutEflow,
            event.MET_CellOutCorr_STVF_etx,
            event.MET_CellOutCorr_STVF_ety,
            event.MET_CellOutCorr_STVF_sumet)

        MET = self.met_utility.getMissingET(METUtil.RefFinal)

        if self.verbose:
            print "Recalculated MET: %.3f (original: %.3f)" % (
                    MET.et(), event.MET_RefFinal_STVF_et)

        # update the MET with the shifted value
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

"""
Adapted from the Example.C in MissingETUtility/macros
"""
# stdlib Python imports
import sys
import math
from math import sin, sqrt, pow

from atlastools import utils
from atlastools import datasets

# ROOT imports
import ROOT
from ROOT import TH1D, TFile

# ATLAS tools imports
from externaltools import MissingETUtility
from externaltools import MuonMomentumCorrections
from externaltools import JetUncertainties
from externaltools import JetResolution
from externaltools import egammaAnalysisUtils

# MissingETUtility
from ROOT import METUtility
from ROOT import METUtil
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


class Systematics(EventFilter):

    # electrons
    EES_UP = METUtil.EESUp
    EES_DOWN = METUtil.EESDown
    EER_UP = METUtil.EERUp
    EER_DOWN = METUtil.EERDown

    # photons
    PES_UP = METUtil.PESUp
    PES_DOWN = METUtil.PESDown
    PER_UP = METUtil.PERUp
    PER_DOWN = METUtil.PERDown

    # taus
    TES_UP = METUtil.TESUp
    TES_DOWN = METUtil.TESDown
    TER_UP = METUtil.TERUp
    TER_DOWN = METUtil.TERDown

    # jets
    JES_UP = METUtil.JESUp
    JES_DOWN = METUtil.JESDown
    JER_UP = METUtil.JERUp
    JER_DOWN = METUtil.JERDown

    # muons
    MERID_UP = METUtil.MERIDUp
    MERID_DOWN = METUtil.MERIDDown
    MERMS_UP = METUtil.MERMSUp
    MERMS_DOWN = METUtil.MERMSDown
    MES_UP = METUtil.MESUp
    MES_DOWN = METUtil.MESDown

    # cells ?
    CES_UP = METUtil.CESUp
    CES_DOWN = METUtil.CESDown
    CER_UP = METUtil.CERUp
    CER_DOWN = METUtil.CERDown

    # tracks
    TRKES_UP = METUtil.TrkESUp
    TRKES_DOWN = METUtil.TrkESDown
    TRKER_UP = METUtil.TrkERUp
    TRKER_DOWN = METUtil.TrkERDown

    # clusters
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
            systematic,
            datatype,
            year,
            verbose=False,
            stream=None,
            **kwargs):

        super(Systematics, self).__init__(**kwargs)
        self.systematic = systematic
        self.datatype = datatype
        self.year = year
        self.verbose = verbose
        if stream is None:
            self.stream = sys.stdout
        else:
            self.stream = stream

        # Initialise your METUtility object
        self.testUtil = METUtility()

        # *** Demonstration of the configuration here  ***
        # *** All values that are set are the defaults ***

        # Turn on (off) the relevant MET terms
        # Standard MET_RefFinal has:
        # RefEle, RefGamma, RefTau, RefJet, SoftJets, RefMuon, MuonTotal, (CellOut), CellOut_Eflow
        self.testUtil.defineMissingET(True, True, True, True, True, True, True, False, True)

        # SUSY group, MET_Simplified20 has
        # RefEle, (RefGamma), (RefTau), RefJet, (SoftJets), (RefMuon), MuonTotal, CellOut, (CellOut_Eflow)
        # self.testUtil.defineMissingET(True, False, True, True, False, False, True, True, False)

        # The threshold below which jets enter the SoftJets term (JES is not applied)
        self.testUtil.setSoftJetCut(20e3)

        # Whether to use MUID muons (otherwise STACO).
        self.testUtil.setIsMuid(False)

        # Whether METUtility should scream at you over every little thing
        self.testUtil.setVerbosity(False)

        self.systUtil = METUtility()

        # Some other options are available, but not recommended/deprecated

        # *** Set up the uncertainty tools ***
        # Tag assumed: JetUncertainties-00-05-09-02
        self.jesTool = MultijetJESUncertaintyProvider(
            "MultijetJES_Preliminary.config",
            "InsituJES2011_AllNuisanceParameters.config",
            "AntiKt4LCTopoJets","MC11c")

        # Tag assumed: JetResolution-01-00-00
        self.jerTool = JERProvider(
            "AntiKt4LCTopoJES", "Truth",
            JetResolution.get_resource('JERProviderPlots.root'))
        self.jerTool.init()

        # Tag assumed: egammaAnalysisUtils-00-02-76
        self.egammaTool = eg2011.EnergyRescaler()
        self.egammaTool.useDefaultCalibConstants("2011")

        # Tag assumed: MuonMomentumCorrections-00-05-03
        self.muonTool = MuonSmear.SmearingClass("Data11","staco","q_pT","Rel17")
        self.muonTool.UseScale(1)
        self.muonTool.UseImprovedCombine()

        # No tag yet, testing code
        self.tesTool = TESUncertaintyProvider()

    def passes(self, event):

        # ResoSoftTerms uses gRandom for smearing. Set the seed here however you like.
        if self.datatype = datasets.DATA:
            ROOT.gRandom.SetSeed(int(event.RunNumber * event.EventNumber))
        else:
            ROOT.gRandom.SetSeed(int(event.mc_channel_number * event.EventNumber))

        # Check for a good primary vertex
        # This is needed for jet and soft term systematics
        goodPV = False
        nvtxsoftmet = 0
        nvtxjets = 0

        if event.vertices:
            # Most D3PDs contain the vx_type branch, but some don't.
            # Those which don't are most likely skimmed, and require at least 1 primary vertex for all events.
            # If your D3PD is skimmed in this way, then the goodPV (nTracks and z) check should be applied
            # to the first vertex in the collection.
            # Otherwise, you should ensure that the vx_type branch is available.
            for vertex in event.vertices:
                if vertex.type == 1 and vertex.nTracks > 2 and abs(vertex.z) < 200:
                    goodPV = True
            if goodPV:
                for vertex in event.vertices:
                    if vertex.nTracks > 2:
                        nvtxsoftmet += 1
                    if vertex.nTracks > 1:
                        nvtxjets += 1

        # See this method for an example of how to pass object systematics
        # to METUtility such that they can be propagated to the MET.
        self.run_systematics()

        # This demonstration is for doing smearing and systematics
        self.systUtil.reset()
        self.systUtil.setJetParameters(
                event.jet_pt,
                event.jet_eta,
                event.jet_phi,
                event.jet_E,
                event.jet_MET_wet,
                event.jet_MET_wpx,
                event.jet_MET_wpy,
                event.jet_MET_statusWord)

        self.systUtil.setOriJetParameters(event.jet_pt)

        # Putting in smeared and/or scaled objects will cause that to be reflected in MET
        self.systUtil.setElectronParameters(
                event.el_smeared_pt,
                event.el_eta,
                event.el_phi,
                event.el_MET_wet,
                event.el_MET_wpx,
                event.el_MET_wpy,
                event.el_MET_statusWord)

        self.systUtil.setPhotonParameters(
                event.ph_smeared_pt,
                event.ph_eta,
                event.ph_phi,
                event.ph_MET_wet,
                event.ph_MET_wpx,
                event.ph_MET_wpy,
                event.ph_MET_statusWord)

        self.systUtil.setTauParameters(
                event.tau_pt,
                event.tau_eta,
                event.tau_phi,
                event.tau_MET_wet,
                event.tau_MET_wpx,
                event.tau_MET_wpy,
                event.tau_MET_statusWord)

        self.systUtil.setMuonParameters(
                event.mu_smeared_pt,
                event.mu_eta,
                event.mu_phi,
                event.mu_MET_wet,
                event.mu_MET_wpx,
                event.mu_MET_wpy,
                event.mu_MET_statusWord)

        # In this instance there is an overloaded version of setExtraMuonParameters that accepts smeared pTs for spectro
        self.systUtil.setExtraMuonParameters(
                event.mu_smeared_ms_pt,
                event.mu_ms_theta,
                event.mu_ms_phi)

        self.systUtil.setMETTerm(
                METUtil.RefTau,
                event.MET_RefTau_etx,
                event.MET_RefTau_ety,
                event.MET_RefTau_sumet)

        self.systUtil.setMETTerm(
                METUtil.RefMuon,
                event.MET_RefMuon_etx,
                event.MET_RefMuon_ety,
                event.MET_RefMuon_sumet)

        self.systUtil.setMETTerm(
                METUtil.SoftJets,
                event.MET_SoftJets_etx,
                event.MET_SoftJets_ety,
                event.MET_SoftJets_sumet)

        #   self.systUtil.setMETTerm(METUtil.CellOut, MET_CellOut_etx, MET_CellOut_ety, MET_CellOut_sumet)
        self.systUtil.setMETTerm(
                METUtil.CellOutEflow,
                event.MET_CellOut_Eflow_etx,
                event.MET_CellOut_Eflow_ety,
                event.MET_CellOut_Eflow_sumet)

        # These set up the systematic "SoftTerms_ptHard"
        self.systUtil.setNvtx(nvtxsoftmet)
        #  if self.datatype != datasets.DATA:
        #    self.systUtil.setMETTerm(METUtil.Truth, MET_Truth_NonInt_etx, MET_Truth_NonInt_ety, MET_Truth_NonInt_sumet)

        self.systUtil.setObjectEnergyUncertainties(
                METUtil.Jets,
                jesUp,
                jesDown)

        self.systUtil.setObjectResolutionShift(
                METUtil.Jets,
                jerUp,
                jerDown)

        self.systUtil.setObjectEnergyUncertainties(
                METUtil.Electrons,
                eesUp,
                eesDown)

        self.systUtil.setObjectResolutionShift(
                METUtil.Electrons,
                eerUp,
                eerDown)

        self.systUtil.setObjectEnergyUncertainties(
                METUtil.Photons,
                pesUp,
                pesDown)

        self.systUtil.setObjectResolutionShift(
                METUtil.Photons,
                perUp,
                perDown)

        # Muons are complicated, and MET makes use of track, spectro, and combined quantites,
        # so we need all of their resolutions.
        # comboms reflects that it is the combined muon res affected by shifting ms res up and down.
        # comboid is for shifting the id res up and down
        self.systUtil.setObjectResolutionShift(
                METUtil.MuonsComboMS,
                cb_mermsUp,
                cb_mermsDown)

        self.systUtil.setObjectResolutionShift(
                METUtil.MuonsComboID,
                cb_meridUp,
                cb_meridDown)

        self.systUtil.setObjectResolutionShift(
                METUtil.SpectroMuons,
                mermsUp,
                mermsDown)

        # For now the mes only affects combined
        self.systUtil.setObjectEnergyUncertainties(
                METUtil.Muons,
                mesUp,
                mesDown)

        # Taus are just in as an example
        self.systUtil.setObjectEnergyUncertainties(
                METUtil.Taus,
                tesUp,
                tesDown)

        # Fill histograms
        met_RefFinal = self.systUtil.getMissingET(METUtil.RefFinal)
        self.h_RefFinal.Fill(met_RefFinal.et()/1000)

        # Jet systematics
        met_RefFinal_JESUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.JESUp)
        self.h_RefFinal_JESUp.Fill(met_RefFinal_JESUp.et()/1000)

        met_RefFinal_JESDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.JESDown)
        self.h_RefFinal_JESDown.Fill(met_RefFinal_JESDown.et()/1000)

        met_RefFinal_JERUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.JERUp)
        self.h_RefFinal_JERUp.Fill(met_RefFinal_JERUp.et()/1000)

        # Electron systematics
        met_RefFinal_EESUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.EESUp)
        self.h_RefFinal_EESUp.Fill(met_RefFinal_EESUp.et()/1000)

        met_RefFinal_EESDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.EESDown)
        self.h_RefFinal_EESDown.Fill(met_RefFinal_EESDown.et()/1000)

        met_RefFinal_EERUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.EERUp)
        self.h_RefFinal_EERUp.Fill(met_RefFinal_EERUp.et()/1000)

        met_RefFinal_EERDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.EERDown)
        self.h_RefFinal_EERDown.Fill(met_RefFinal_EERDown.et()/1000)

        # Photon systematics
        met_RefFinal_PESUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.PESUp)
        self.h_RefFinal_PESUp.Fill(met_RefFinal_PESUp.et()/1000)

        met_RefFinal_PESDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.PESDown)
        self.h_RefFinal_PESDown.Fill(met_RefFinal_PESDown.et()/1000)

        met_RefFinal_PERUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.PERUp)
        self.h_RefFinal_PERUp.Fill(met_RefFinal_PERUp.et()/1000)

        met_RefFinal_PERDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.PERDown)
        self.h_RefFinal_PERDown.Fill(met_RefFinal_PERDown.et()/1000)

        # Muon systematics
        met_RefFinal_MESUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.MESUp)
        self.h_RefFinal_MESUp.Fill(met_RefFinal_MESUp.et()/1000)

        met_RefFinal_MESDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.MESDown)
        self.h_RefFinal_MESDown.Fill(met_RefFinal_MESDown.et()/1000)

        met_RefFinal_MERIDUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.MERIDUp)
        self.h_RefFinal_MERIDUp.Fill(met_RefFinal_MERIDUp.et()/1000)

        met_RefFinal_MERIDDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.MERIDDown)
        self.h_RefFinal_MERIDDown.Fill(met_RefFinal_MERIDDown.et()/1000)

        met_RefFinal_MERMSUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.MERMSUp)
        self.h_RefFinal_MERMSUp.Fill(met_RefFinal_MERMSUp.et()/1000)

        met_RefFinal_MERMSDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.MERMSDown)
        self.h_RefFinal_MERMSDown.Fill(met_RefFinal_MERMSDown.et()/1000)

        # Soft terms systematics
        met_RefFinal_ScaleSoftTermsUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ScaleSoftTermsUp)
        self.h_RefFinal_ScaleSoftTermsUp.Fill(met_RefFinal_ScaleSoftTermsUp.et()/1000)

        met_RefFinal_ScaleSoftTermsDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ScaleSoftTermsDown)
        self.h_RefFinal_ScaleSoftTermsDown.Fill(met_RefFinal_ScaleSoftTermsDown.et()/1000)

        met_RefFinal_ResoSoftTermsUp = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ResoSoftTermsUp)
        self.h_RefFinal_ResoSoftTermsUp.Fill(met_RefFinal_ResoSoftTermsUp.et()/1000)

        met_RefFinal_ResoSoftTermsDown = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ResoSoftTermsDown)
        self.h_RefFinal_ResoSoftTermsDown.Fill(met_RefFinal_ResoSoftTermsDown.et()/1000)

        # Alternate soft terms systematics
        met_RefFinal_ScaleSoftTermsUp_ptHard = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ScaleSoftTermsUp_ptHard)
        self.h_RefFinal_ScaleSoftTermsUp_ptHard.Fill(met_RefFinal_ScaleSoftTermsUp_ptHard.et()/1000)

        met_RefFinal_ScaleSoftTermsDown_ptHard = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ScaleSoftTermsDown_ptHard)
        self.h_RefFinal_ScaleSoftTermsDown_ptHard.Fill(met_RefFinal_ScaleSoftTermsDown_ptHard.et()/1000)

        met_RefFinal_ResoSoftTermsUp_ptHard = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ResoSoftTermsUp_ptHard)
        self.h_RefFinal_ResoSoftTermsUp_ptHard.Fill(met_RefFinal_ResoSoftTermsUp_ptHard.et()/1000)

        met_RefFinal_ResoSoftTermsUpDown_ptHard = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ResoSoftTermsUpDown_ptHard)
        self.h_RefFinal_ResoSoftTermsUpDown_ptHard.Fill(met_RefFinal_ResoSoftTermsUpDown_ptHard.et()/1000)

        met_RefFinal_ResoSoftTermsDownUp_ptHard = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ResoSoftTermsDownUp_ptHard)
        self.h_RefFinal_ResoSoftTermsDownUp_ptHard.Fill(met_RefFinal_ResoSoftTermsDownUp_ptHard.et()/1000)

        met_RefFinal_ResoSoftTermsDown_ptHard = self.systUtil.getMissingET(METUtil.RefFinal,METUtil.ResoSoftTermsDown_ptHard)
        self.h_RefFinal_ResoSoftTermsDown_ptHard.Fill(met_RefFinal_ResoSoftTermsDown_ptHard.et()/1000)
        return True

    def DemoMetSystematics(self, event):
        #######################################
        # Demonstrates how to set up the METUtility with object momenta
        # such that MET_RefFinal can be rebuilt matching the values in D3PD.
        #
        # *** *** *** *** *** DISCLAIMER *** *** *** *** ***
        #
        # These examples of uncertainty-setting are meant to
        # demonstrate how the uncertainties should be passed
        # to METUtility.  Recommendations on just which ones
        # you are meant to be using come from the CP groups.


        """
        JET SYSTEMATICS
        """
        # First, we get the jet energy scale uncertainties and
        # systematic variation in the jet resolutions
        jesUp = ROOT.vector('float')()
        jesDown = ROOT.vector('float')()
        jerUp = ROOT.vector('float')()
        jerDown = ROOT.vector('float')()

        # Note on use of ROOT random number generators:
        # TRandom and TRandom2 have many documented deficiencies.
        # TRandom3 is generally considered safely usable.
        # Also note that ROOT's gRandom calls TRandom3.
        jetRandom = ROOT.TRandom3()

        # jet_ is jet_AntiKt4LCTopo in tau D3PDs
        for ijet, jet in enumerate(event.jets):

            jesShiftUp = 0.0
            jesShiftDown = 0.0
            jerShift = 1.0
            # Safest to assume nothing about the uncertainties on soft jets.
            # These will go into SoftJets anyhow, and so the JES systematics
            # aren't used.

            if jet.pt > 20e3 and jet.pt < 7000e3 and abs(jet.eta) < 4.5:

                # delta R cut needed to apply close-by jets uncertainty
                drmin=9999
                pi = math.pi
                for jjet, otherjet in enumerate(event.jets):
                    if jjet.emscale_pt > 7000:
                        if ijet != jjet:
                            dr = utils.dR(jet.eta, jet.phi,
                                          otherjet.eta,
                                          otherjet.phi)
                            if dr < drmin:
                                drmin = dr

                # The bool is the "isPos" argument
                jesShiftUp = self.jesTool.getRelUncert(jet.pt,
                                 jet.eta,drmin,
                                 True, nvtxjets, averageIntPerXing)
                jesShiftDown = -1*self.jesTool.getRelUncert(jet.pt,
                                  jet.eta,drmin,
                                  False, nvtxjets, averageIntPerXing)

            jesUp.push_back(jesShiftUp)
            jesDown.push_back(jesShiftDown)

            # Allowable range is > 10 GeV, but anything below 20 enters SoftJets
            if jet.pt > 20e3 and jet.pt < 5000e3:
                pt = jet.pt
                eta = jet.eta
                if abs(eta) > 4.5:
                    eta = 4.49 if eta > 0 else -4.49

                S = self.jerTool.getRelResolutionMC(pt/1e3,eta)
                U = self.jerTool.getResolutionUncert(pt/1e3,eta)
                smearingFactorSyst = sqrt(pow(S+U,2)-pow(S,2))

                # You can set the seed however you like, but if reproducibility
                # is required, setting it to something like object phi ensures
                # a good mix of randomness and reproducibility.
                jetRandom.SetSeed(int(1.e5 * jet.phi))
                jerShift = jetRandom.Gaus(0, smearingFactorSyst)

            jerUp.push_back(jerShift)
            jerDown.push_back(-1 * jerShift); # Usually not used, see below.

            ###################################
            # Note: The JERDown shift is essentially meaningless.
            # If one is smearing central values, then there is an alternate
            # definition, i.e. from r16:
            #
            # S = self.jerTool.getRelResolutionData(pt/1e3,eta)
            # SMC = self.jerTool.getRelResolutionMC(pt/1e3,eta)
            # U = self.jerTool.getResolutionUncert(pt/1e3,eta)
            # smearingFactorMC = sqrt( S*S - SMC*SMC )
            # smearingFactorSystUp = sqrt( (S+U)*(S+U) - SMC*SMC )
            # smearingFactorSystDown = (S-U > SMC) ? sqrt( (S+U)*(S+U) - SMC*SMC ) : 0
            #
            # float jerShift = jetRandom.Gaus(1,smearingFactorMC)
            # float jerShiftUp = jetRandom.Gaus(1,smearingFactorSystUp)/jerShift
            # float jerShiftDown = jetRandom.Gaus(1,smearingFactorSystDown)/jerShift
            #
            # jet_smeared_pt = pt*jerShift
            # jerUp.push_back(jerShiftUp-1)
            # jerDown.push_back(jerShiftDown-1)
            #
            # This means that we smear the MC jets to match the resolution in data
            # for central values, or the resolution +/- uncertainty.
            # The standard practice is only to use res + uncertainty.
            #
            ###################################

        del jetRandom

        """
        ELECTRON SYSTEMATICS
        """
        # Here we get the electron energy scale and resolution systematics
        eesUp = ROOT.vector('float')()
        eesDown = ROOT.vector('float')()
        eerUp = ROOT.vector('float')()
        eerDown = ROOT.vector('float')()
        el_smeared_pt = ROOT.vector('float')()

        for iel, el in enumerate(event.electrons):

            self.egammaTool.SetRandomSeed(int(1e5*abs(el.phi)))

            # Smear to match the data resolution, or by systematic variations
            smear = self.egammaTool.getSmearingCorrectionMeV(el.cl_eta, el.E, 0, True)
            smearUp = self.egammaTool.getSmearingCorrectionMeV(el.cl_eta, el.E, 2, True)
            smearDown = self.egammaTool.getSmearingCorrectionMeV(el.cl_eta, el.E, 1, True)

            el_smeared_pt.push_back(smear * el.pt)
            eerUp.push_back((smearUp - smear) / smear)
            eerDown.push_back((smearDown - smear) / smear)

            # Correct the measured energies in data, and scale by systematic variations
            correction = 1.
            if self.datatype == datasets.DATA:
                correction = self.egammaTool.applyEnergyCorrectionMeV(
                        el.cl_eta,
                        el.cl_phi,
                        el.E,
                        el.cl_pt,
                        0,"ELECTRON") / el.E

            el_smeared_pt.at(iel) *= correction
            energyUp = self.egammaTool.applyEnergyCorrectionMeV(
                    el.cl_eta,
                    el.cl_phi,
                    el.E,
                    el.cl_pt,
                    2,"ELECTRON") / (correction * el.E) - 1
            energyDown = self.egammaTool.applyEnergyCorrectionMeV(
                    el.cl_eta,
                    el.cl_phi,
                    el.E,
                    el.cl_pt,
                    1,"ELECTRON") / (correction * el.E) - 1

            eesUp.push_back(energyUp)
            eesDown.push_back(energyDown)

        """
        PHOTON SYSTEMATICS
        """
        # Now we get the same for photons
        pesUp = ROOT.vector('float')()
        pesDown = ROOT.vector('float')()
        perUp = ROOT.vector('float')()
        perDown = ROOT.vector('float')()
        ph_smeared_pt = ROOT.vector('float')()

        for iph, ph in enumerate(event.photons):

            self.egammaTool.SetRandomSeed(int(1e5*abs(ph.phi)))

            # Smear to match the data resolution, or by systematic variations
            smear = self.egammaTool.getSmearingCorrectionMeV(ph.cl_eta, ph.E, 0, True)
            smearUp = self.egammaTool.getSmearingCorrectionMeV(ph.cl_eta, ph.E, 2, True)
            smearDown = self.egammaTool.getSmearingCorrectionMeV(ph.cl_eta, ph.E, 1, True)

            ph_smeared_pt.push_back(smear * ph.pt)
            perUp.push_back((smearUp - smear) / smear)
            perDown.push_back((smearDown - smear) / smear)

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

            ph_smeared_pt.at(iph) *= correction
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

            pesUp.push_back(energyUp)
            pesDown.push_back(energyDown)

        """
        MUON SYSTEMATICS
        """
        # And now the same for muons. We need resolution shifts for ID and MS,
        # and different treatment for the MS four-vector (for standalone muons).
        mu_smeared_pt = ROOT.vector('float')()
        mu_smeared_ms_pt = ROOT.vector('float')()

        cb_meridUp = ROOT.vector('float')()
        cb_meridDown = ROOT.vector('float')()
        cb_mermsUp = ROOT.vector('float')()
        cb_mermsDown = ROOT.vector('float')()
        mermsUp = ROOT.vector('float')()
        mermsDown = ROOT.vector('float')()

        mesUp = ROOT.vector('float')()
        mesDown = ROOT.vector('float')()

        for mu in event.muons:

            ptcb = mu.pt
            ptid = (mu.id_qoverp_exPV != 0.) ? abs(sin(mu.id_theta_exPV)/mu.id_qoverp_exPV) : 0.
            ptms = (mu.ms_qoverp != 0.) ? abs(sin(mu.ms_theta)/mu.ms_qoverp) : 0.
            self.muonTool.SetSeed(int(1.e+5*fabs(mu.phi)))
            etaMu = mu.eta
            charge = mu.charge
            self.muonTool.Event(ptms, ptid, ptcb, etaMu, charge)

            smearedCombinedPt = self.muonTool.pTCB()
            if not mu.isCombinedMuon:
                smearedCombinedPt = self.muonTool.pTMS() + self.muonTool.pTID()

            smearedMSPt = self.muonTool.pTMS()

            mu_smeared_ms_pt.push_back(smearedMSPt)
            mu_smeared_pt.push_back(smearedCombinedPt)

            smearedpTMS = 0.1
            smearedpTID = 0.1
            smearedpTCB = 0.1

            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSLOW")
            smearedpTMS = ptMS_smeared/smearedMSPt - 1.0
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            mermsDown.push_back(smearedpTMS)
            cb_mermsDown.push_back(smearedpTCB)
            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSUP")
            smearedpTMS = ptMS_smeared/smearedMSPt - 1.0
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            mermsUp.push_back(smearedpTMS)
            cb_mermsUp.push_back(smearedpTCB)
            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDUP")
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            cb_meridUp.push_back(smearedpTCB)
            self.muonTool.PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDLOW")
            smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0
            cb_meridDown.push_back(smearedpTCB)

            detRegion = self.muonTool.DetRegion()
            if detRegion == -1:
                detRegion = 3
            scalesyst = self.muonTool.getScaleSyst_CB().at(detRegion)
            mesUp.push_back(scalesyst)
            mesDown.push_back(-scalesyst)

        """
        TAU SYSTEMATICS
        """
        tesUp = ROOT.vector('float')()
        tesDown = ROOT.vector('float')()

        # And for taus (this is test code, do not use without tau group approval)
        for tau in event.taus:
            pt = tau.pt
            eta = tau.eta
            nProng = tau.nProng
            double uncert = self.tesTool.GetTESUncertainty(pt / 1e3, eta, nProng)
            if uncert < 0:
                uncert = 0
            tesUp.push_back(uncert)
            tesDown.push_back(-1 * uncert)



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


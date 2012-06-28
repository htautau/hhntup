"""
Adapted from the Example.C in MissingETUtility/macros
"""
# stdlib Python imports
import sys

# ROOT imports
import ROOT
from ROOT import TH1D, TFile

# ATLAS tools imports
# MissingETUtility
from ROOT import METUtility
# JetUncertainties
from ROOT import MultijetJESUncertaintyProvider
# JetResolution
from ROOT import JERProvider
# egammaAnalysisUtils
from ROOT import EnergyRescaler
# MuonMomentumCorrections
from ROOT import SmearingClass
# MissingETUtility
from ROOT import TESUncertaintyProvider


class Systematic(EventFilter):

    def __init__(self, datatype, year, verbose=False, stream=None, **kwargs):

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
            "../../Jet/JetResolution/share/JERProviderPlots.root")
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

        # Prepare an output file
        self.outfile = TFile("metutil_example.root","RECREATE")
        # Set up a bunch of histograms
        h_RefFinal = TH1D("h_RefFinal", "RefFinal",50,0,1000)
        # Jet systematics
        h_RefFinal_JESUp = TH1D("h_RefFinal_JESUp", "RefFinal JESUp",50,0,1000)
        h_RefFinal_JESDown = TH1D("h_RefFinal_JESDown", "RefFinal JESDown",50,0,1000)
        h_RefFinal_JERUp = TH1D("h_RefFinal_JERUp", "RefFinal JERUp",50,0,1000)
        # Electron systematics
        h_RefFinal_EESUp = TH1D("h_RefFinal_EESUp", "RefFinal EESUp",50,0,1000)
        h_RefFinal_EESDown = TH1D("h_RefFinal_EESDown",
        "RefFinal EESDown",50,0,1000)
        h_RefFinal_EERUp = TH1D("h_RefFinal_EERUp",
          "RefFinal EERUp",50,0,1000)
        h_RefFinal_EERDown = TH1D("h_RefFinal_EERDown",
        "RefFinal EERDown",50,0,1000)
        # Photon systematics
        h_RefFinal_PESUp = TH1D("h_RefFinal_PESUp",
          "RefFinal PESUp",50,0,1000)
        h_RefFinal_PESDown = TH1D("h_RefFinal_PESDown",
        "RefFinal PESDown",50,0,1000)
        h_RefFinal_PERUp = TH1D("h_RefFinal_PERUp",
          "RefFinal PERUp",50,0,1000)
        h_RefFinal_PERDown = TH1D("h_RefFinal_PERDown",
        "RefFinal PERDown",50,0,1000)
        # Muon systematics
        h_RefFinal_MESUp = TH1D("h_RefFinal_MESUp",
          "RefFinal MESUp",50,0,1000)
        h_RefFinal_MESDown = TH1D("h_RefFinal_MESDown",
        "RefFinal MESDown",50,0,1000)
        h_RefFinal_MERIDUp = TH1D("h_RefFinal_MERIDUp",
        "RefFinal MERIDUp",50,0,1000)
        h_RefFinal_MERIDDown = TH1D("h_RefFinal_MERIDDown",
          "RefFinal MERIDDown",50,0,1000)
        h_RefFinal_MERMSUp = TH1D("h_RefFinal_MERMSUp",
        "RefFinal MERMSUp",50,0,1000)
        h_RefFinal_MERMSDown = TH1D("h_RefFinal_MERMSDown",
          "RefFinal MERMSDown",50,0,1000)
        # Soft terms systematics
        h_RefFinal_ScaleSoftTermsUp = TH1D("h_RefFinal_ScaleSoftTermsUp",
             "RefFinal ScaleSoftTermsUp",50,0,1000)
        h_RefFinal_ScaleSoftTermsDown = TH1D("h_RefFinal_ScaleSoftTermsDown",
               "RefFinal ScaleSoftTermsDown",50,0,1000)
        h_RefFinal_ResoSoftTermsUp = TH1D("h_RefFinal_ResoSoftTermsUp",
            "RefFinal ResoSoftTermsUp",50,0,1000)
        h_RefFinal_ResoSoftTermsDown = TH1D("h_RefFinal_ResoSoftTermsDown",
              "RefFinal ResoSoftTermsDown",50,0,1000)
        # Alternate soft terms systematics
        h_RefFinal_ScaleSoftTermsUp_ptHard = TH1D("h_RefFinal_ScaleSoftTermsUp_ptHard",
                "RefFinal ScaleSoftTermsUp_ptHard",50,0,1000)
        h_RefFinal_ScaleSoftTermsDown_ptHard = TH1D("h_RefFinal_ScaleSoftTermsDown_ptHard",
                  "RefFinal ScaleSoftTermsDown_ptHard",50,0,1000)
        h_RefFinal_ResoSoftTermsUp_ptHard = TH1D("h_RefFinal_ResoSoftTermsUp_ptHard",
                   "RefFinal ResoSoftTermsUp_ptHard",50,0,1000)
        h_RefFinal_ResoSoftTermsUpDown_ptHard = TH1D("h_RefFinal_ResoSoftTermsUpDown_ptHard",
                   "RefFinal ResoSoftTermsUpDown_ptHard",50,0,1000)
        h_RefFinal_ResoSoftTermsDownUp_ptHard = TH1D("h_RefFinal_ResoSoftTermsDownUp_ptHard",
                   "RefFinal ResoSoftTermsDownUp_ptHard",50,0,1000)
        h_RefFinal_ResoSoftTermsDown_ptHard = TH1D("h_RefFinal_ResoSoftTermsDown_ptHard",
                 "RefFinal ResoSoftTermsDown_ptHard",50,0,1000)


    def passes(self, event):

        # See this method for an example of setting up the METUtility
        # such that it can rebuild the D3PD nominal MET values.
        self.CheckMetRebuilding()

        # See this method for an example of how to pass object systematics
        # to METUtility such that they can be propagated to the MET.
        self.DemoMetSystematics()
        return True

    def CheckMetRebuilding(self, event):

        # Demonstrates how to set up the METUtility with object momenta such that
        # MET_RefFinal can be rebuilt matching the values in D3PD.

        # Start with a clean METUtility
        self.testUtil.reset()

        # For these, just the kinematics need to be set
        self.testUtil.setElectronParameters(el_pt, el_eta, el_phi,
                        el_MET_wet, el_MET_wpx, el_MET_wpy,
                        el_MET_statusWord)

        self.testUtil.setPhotonParameters(ph_pt, ph_eta, ph_phi,
                      ph_MET_wet, ph_MET_wpx, ph_MET_wpy,
                      ph_MET_statusWord)

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
        self.testUtil.setJetParameters(jet_pt, jet_eta, jet_phi, jet_E,
                       jet_MET_wet, jet_MET_wpx, jet_MET_wpy,
                       jet_MET_statusWord)
        self.testUtil.setOriJetParameters(jet_pt)

        # Muons may be ID, combined, or standalone. For the latter especially,
        # we need to set the MS four-momenta because they are treated differently.
        self.testUtil.setMuonParameters(mu_pt, mu_eta, mu_phi,
                      mu_MET_wet, mu_MET_wpx, mu_MET_wpy,
                      mu_MET_statusWord)
        self.testUtil.setExtraMuonParameters(mu_ms_qoverp, mu_ms_theta, mu_ms_phi, mu_charge)
        # An alternative version of this method is demonstrated below, and takes pT, eta, phi instead.
        # This is more convenient when one needs to smear the pT, for example.

        # When the terms are not rebuilt from the objects, due to incomplete info,
        # then one needs to set the term directly from the stored value in the D3PD.
        # This might also be done if you aren't interested in the possible variations
        # of that term. E.g. if you only care about photons, no need to rebuild muons.

        # self.testUtil.setMETTerm(METUtil::RefJet, MET_RefJet_etx, MET_RefJet_ety, MET_RefJet_sumet)
        self.testUtil.setMETTerm(METUtil::SoftJets, MET_SoftJets_etx, MET_SoftJets_ety, MET_SoftJets_sumet)
        # self.testUtil.setMETTerm(METUtil::RefEle, MET_RefEle_etx, MET_RefEle_ety, MET_RefEle_sumet)
        # self.testUtil.setMETTerm(METUtil::RefGamma, MET_RefGamma_etx, MET_RefGamma_ety, MET_RefGamma_sumet)

        # *** Note the difference in naming -- there is normally no MET_MuonTotal term.
        #     It's usually Muon_Total, Muon_Total_Muid, something like that.
        #     MET_RefFinal in particular uses MET_MuonBoy.
        # self.testUtil.setMETTerm(METUtil::MuonTotal, MET_MuonBoy_etx, MET_MuonBoy_ety, MET_MuonBoy_sumet)

        # *** Note that RefMuon is not rebuilt from muons -- it is a calorimeter term.
        self.testUtil.setMETTerm(METUtil::RefMuon, MET_RefMuon_etx, MET_RefMuon_ety, MET_RefMuon_sumet)
        self.testUtil.setMETTerm(METUtil::RefTau, MET_RefTau_etx, MET_RefTau_ety, MET_RefTau_sumet)
        # self.testUtil.setMETTerm(METUtil::CellOut, MET_CellOut_etx, MET_CellOut_ety, MET_CellOut_sumet)
        self.testUtil.setMETTerm(METUtil::CellOutEflow, MET_CellOut_Eflow_etx, MET_CellOut_Eflow_ety, MET_CellOut_Eflow_sumet)

        # This is the simple check, where you compare the terms manually against what's in the D3PD.
        # Note: every call to getMissingET recomputes the terms, so if you need to get more than one
        # value, e.g. etx, ety, et, sumET, it's more efficient to get the METObject.
        # Usually, comparing etx and/or ety is more informative, because et could be right if
        # etx and ety were flipped, for example. They also add linearly, which et doesn't.

        METObject RefEle_util = self.testUtil.getMissingET(METUtil::RefEle)
        METObject RefGamma_util = self.testUtil.getMissingET(METUtil::RefGamma)
        METObject RefTau_util = self.testUtil.getMissingET(METUtil::RefTau)
        METObject RefMuon_util = self.testUtil.getMissingET(METUtil::RefMuon)
        METObject RefJet_util = self.testUtil.getMissingET(METUtil::RefJet)
        METObject SoftJets_util = self.testUtil.getMissingET(METUtil::SoftJets)
        #      METObject refCellOut_util = self.testUtil.getMissingET(METUtil::CellOut)
        METObject CellOutEflow_util = self.testUtil.getMissingET(METUtil::CellOutEflow)
        METObject MuonTotal_util = self.testUtil.getMissingET(METUtil::MuonTotal)
        METObject RefFinal_util = self.testUtil.getMissingET(METUtil::RefFinal)

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
        check_refFinal = self.testUtil.checkConsistency(refFinal_test,METUtil::RefFinal)
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
        # By default, METObject::sig() returns etx() / ( 0.5*sqrt(sumet()) )
        #
        # There is also the possibility of returning a more sophisticated estimator for
        # the significance, activated by calling METUtility::doSignificance() in setup.
        # This is still under development, and requires all object resolutions to be set.

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
        for jet in event.jets:
            jesShiftUp = 0.0
            jesShiftDown = 0.0
            jerShift = 1.0
            # Safest to assume nothing about the uncertainties on soft jets.
            # These will go into SoftJets anyhow, and so the JES systematics
            # aren't used.

            if(jet_pt.at(iJet) > 20e3
               && jet_pt.at(iJet) < 7000e3
               && fabs(jet_eta.at(iJet)) < 4.5){

              # delta R cut needed to apply close-by jets uncertainty
              float drmin=9999
              double pi = TMath::Pi()
              if( jet_pt.at(iJet)>20000) {
            for (int ii = 0; ii < jet_n; ii++ ) {
              if(jet_emscale_pt.at(ii)>7000) {
                if(iJet!=ii) {
                  double deta = jet_eta.at(iJet) - jet_eta.at(ii)
                  double dphi = fabs(fmod((jet_phi.at(iJet)
                               - jet_phi.at(ii))+3*pi,2*pi)-pi)
                  double dr = sqrt(deta*deta+dphi*dphi)
                  if(dr<drmin) drmin=dr
                }
              }
            }

        # The bool is the "isPos" argument
        jesShiftUp = self.jesTool.getRelUncert(jet_pt.at(iJet),
                         jet_eta.at(iJet),drmin,
                         True, nvtxjets, averageIntPerXing)
        jesShiftDown = -1*self.jesTool.getRelUncert(jet_pt.at(iJet),
                          jet_eta.at(iJet),drmin,
                          False, nvtxjets, averageIntPerXing)
        }
        jesUp.push_back(jesShiftUp)
        jesDown.push_back(jesShiftDown)

        # Allowable range is > 10 GeV, but anything below 20 enters SoftJets
        if(jet_pt.at(iJet) > 20e3 && jet_pt.at(iJet) < 5000e3){
            double pt = jet_pt.at(iJet)
            double eta = jet_eta.at(iJet)
            if(fabs(eta)>4.5) eta = eta>0 ? 4.49 : -4.49

            double S = self.jerTool.getRelResolutionMC(pt/1e3,eta)
            double U = self.jerTool.getResolutionUncert(pt/1e3,eta)
            double smearingFactorSyst = sqrt(pow(S+U,2)-pow(S,2))

            # You can set the seed however you like, but if reproducibility
            # is required, setting it to something like object phi ensures
            # a good mix of randomness and reproducibility.
            jetRandom.SetSeed(int(1.e5*jet_phi.at(iJet)))
            jerShift = jetRandom.Gaus(0, smearingFactorSyst)
        }

        jerUp.push_back(jerShift)
        jerDown.push_back(-1*jerShift); # Usually not used, see below.

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

        }#end of jet loop

        delete jetRandom

        # Here we get the electron energy scale and resolution systematics
        vector<float> eesUp
        vector<float> eesDown
        vector<float> eerUp
        vector<float> eerDown
        vector<float> *el_smeared_pt = new vector<float>

        for (unsigned int iEl = 0; iEl < el_pt.size(); ++iEl) {

        self.egammaTool.SetRandomSeed(int(1e5*fabs(el_phi.at(iEl))))

        # Smear to match the data resolution, or by systematic variations
        float smear = self.egammaTool.getSmearingCorrectionMeV(el_cl_eta.at(iEl), el_E.at(iEl), 0, True)
        float smearUp = self.egammaTool.getSmearingCorrectionMeV(el_cl_eta.at(iEl), el_E.at(iEl), 2, True)
        float smearDown = self.egammaTool.getSmearingCorrectionMeV(el_cl_eta.at(iEl), el_E.at(iEl), 1, True)

        el_smeared_pt.push_back(smear*el_pt.at(iEl))
        eerUp.push_back((smearUp - smear)/smear)
        eerDown.push_back((smearDown - smear)/smear)

        # Correct the measured energies in data, and scale by systematic variations
        float correction = 1.
        if(isData)
          correction = self.egammaTool.applyEnergyCorrectionMeV(el_cl_eta.at(iEl),el_cl_phi.at(iEl),
                                el_E.at(iEl),el_cl_pt.at(iEl),0,"ELECTRON") / el_E.at(iEl)
        el_smeared_pt.at(iEl)*= correction
        double energyUp = self.egammaTool.applyEnergyCorrectionMeV(el_cl_eta.at(iEl),el_cl_phi.at(iEl),
                                   el_E.at(iEl),el_cl_pt.at(iEl),2,"ELECTRON") / (correction*el_E.at(iEl)) - 1
        double energyDown = self.egammaTool.applyEnergyCorrectionMeV(el_cl_eta.at(iEl),el_cl_phi.at(iEl),
                                     el_E.at(iEl),el_cl_pt.at(iEl),1,"ELECTRON") / (correction*el_E.at(iEl)) - 1

        eesUp.push_back(energyUp)
        eesDown.push_back(energyDown)
        } #end of electron loop


        # Now we get the same for photons
        vector<float> pesUp
        vector<float> pesDown
        vector<float> perUp
        vector<float> perDown
        vector<float> *ph_smeared_pt = new vector<float>

        for (unsigned int iPh = 0; iPh < ph_pt.size(); ++iPh) {

        self.egammaTool.SetRandomSeed(int(1.e+5*fabs(ph_phi.at(iPh))))

        # Smear to match the data resolution, or by systematic variations
        float smear = self.egammaTool.getSmearingCorrectionMeV(ph_cl_eta.at(iPh), ph_E.at(iPh), 0, True)
        float smearUp = self.egammaTool.getSmearingCorrectionMeV(ph_cl_eta.at(iPh), ph_E.at(iPh), 2, True)
        float smearDown = self.egammaTool.getSmearingCorrectionMeV(ph_cl_eta.at(iPh), ph_E.at(iPh), 1, True)

        ph_smeared_pt.push_back(smear*ph_pt.at(iPh))
        perUp.push_back((smearUp - smear)/smear)
        perDown.push_back((smearDown - smear)/smear)

        # Correct the measured energies in data, and scale by systematic variations
        # Conversions are treated differently.
        float correction = 1.
        string photontype = ph_isConv.at(iPh) ? "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON"
        if(isData)
          correction = self.egammaTool.applyEnergyCorrectionMeV(ph_cl_eta.at(iPh),ph_cl_phi.at(iPh),ph_E.at(iPh),ph_cl_pt.at(iPh),0,photontype) / ph_E.at(iPh)
        ph_smeared_pt.at(iPh)*= correction
        double energyUp = self.egammaTool.applyEnergyCorrectionMeV(ph_cl_eta.at(iPh),ph_cl_phi.at(iPh),
                                     ph_E.at(iPh),ph_cl_pt.at(iPh),2,photontype) / (correction*ph_E.at(iPh)) - 1
        double energyDown = self.egammaTool.applyEnergyCorrectionMeV(ph_cl_eta.at(iPh),ph_cl_phi.at(iPh),
                                       ph_E.at(iPh),ph_cl_pt.at(iPh),1,photontype) / (correction*ph_E.at(iPh)) - 1

        pesUp.push_back(energyUp)
        pesDown.push_back(energyDown)
        }#end of photon loop

        # And now the same for muons. We need resolution shifts for ID and MS,
        # and different treatment for the MS four-vector (for standalone muons).
        vector<float> *mu_smeared_pt = new vector<float>
        vector<float> *mu_smeared_ms_pt = new vector<float>

        vector<float> cb_meridUp
        vector<float> cb_meridDown
        vector<float> cb_mermsUp
        vector<float> cb_mermsDown
        vector<float> mermsUp
        vector<float> mermsDown

        vector<float> mesUp
        vector<float> mesDown

        for(unsigned int iMu = 0; iMu < mu_pt.size(); ++iMu){

        double ptcb = mu_pt.at(iMu)
        double ptid = (mu_id_qoverp_exPV.at(iMu) != 0.) ? fabs(sin(mu_id_theta_exPV.at(iMu))/mu_id_qoverp_exPV.at(iMu)) : 0.
        double ptms = (mu_ms_qoverp.at(iMu) != 0.) ? fabs(sin(mu_ms_theta.at(iMu))/mu_ms_qoverp.at(iMu)) : 0.
        self.muonTool.SetSeed(int(1.e+5*fabs(mu_phi.at(iMu))))
        double etaMu = mu_eta.at(iMu)
        double charge = mu_charge.at(iMu)
        self.muonTool.Event(ptms, ptid, ptcb, etaMu, charge)

        Float_t smearedCombinedPt = self.muonTool.pTCB()
        if(!mu_isCombinedMuon.at(iMu)) smearedCombinedPt = self.muonTool.pTMS() + self.muonTool.pTID()

        Float_t smearedMSPt = self.muonTool.pTMS()

        mu_smeared_ms_pt.push_back(smearedMSPt)
        mu_smeared_pt.push_back(smearedCombinedPt)

        double ptMS_smeared, ptID_smeared, ptCB_smeared
        float smearedpTMS, smearedpTID, smearedpTCB
        smearedpTMS = 0.1; smearedpTID = 0.1; smearedpTCB = 0.1
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

        int detRegion = self.muonTool.DetRegion()
        if(detRegion==-1) detRegion = 3
        double scalesyst = self.muonTool.getScaleSyst_CB().at(detRegion)
        mesUp.push_back(scalesyst)
        mesDown.push_back(-scalesyst)

        }#end of muon loop

        tesUp = ROOT.vector('float')()
        tesDown = ROOT.vector('float')()

        # And for taus (this is test code, do not use without tau group approval)
        for int iTau=0; iTau<tau_n; iTau++:
            double pt = tau_pt.at(iTau)
            double eta = tau_eta.at(iTau)
            int nProng = tau_nProng.at(iTau)
            double uncert = self.tesTool.GetTESUncertainty(pt/1e3, eta, nProng)

            if uncert < 0:
                uncert = 0
            tesUp.push_back(uncert)
            tesDown.push_back(-1*uncert)

        # This demonstration is for doing smearing and systematics
        self.systUtil.reset()
        self.systUtil.setJetParameters(jet_pt, jet_eta, jet_phi, jet_E,
                       jet_MET_wet, jet_MET_wpx, jet_MET_wpy,
                       jet_MET_statusWord)
        self.systUtil.setOriJetParameters(jet_pt)

        # Putting in smeared and/or scaled objects will cause that to be reflected in MET
        self.systUtil.setElectronParameters(el_smeared_pt, el_eta, el_phi,
                        el_MET_wet, el_MET_wpx, el_MET_wpy,
                        el_MET_statusWord)
        self.systUtil.setPhotonParameters(ph_smeared_pt, ph_eta, ph_phi,
                      ph_MET_wet, ph_MET_wpx, ph_MET_wpy,
                      ph_MET_statusWord)
        self.systUtil.setTauParameters(tau_pt, tau_eta, tau_phi,
                       tau_MET_wet, tau_MET_wpx, tau_MET_wpy,
                       tau_MET_statusWord)
        self.systUtil.setMuonParameters(mu_smeared_pt, mu_eta, mu_phi,
                    mu_MET_wet, mu_MET_wpx, mu_MET_wpy,
                    mu_MET_statusWord)
        # In this instance there is an overloaded version of setExtraMuonParameters that accepts smeared pTs for spectro
        self.systUtil.setExtraMuonParameters(mu_smeared_ms_pt, mu_ms_theta, mu_ms_phi)

        self.systUtil.setMETTerm(METUtil::RefTau, MET_RefTau_etx, MET_RefTau_ety, MET_RefTau_sumet)
        self.systUtil.setMETTerm(METUtil::RefMuon, MET_RefMuon_etx, MET_RefMuon_ety, MET_RefMuon_sumet)
        self.systUtil.setMETTerm(METUtil::SoftJets, MET_SoftJets_etx, MET_SoftJets_ety, MET_SoftJets_sumet)
        #   self.systUtil.setMETTerm(METUtil::CellOut, MET_CellOut_etx, MET_CellOut_ety, MET_CellOut_sumet)
        self.systUtil.setMETTerm(METUtil::CellOutEflow, MET_CellOut_Eflow_etx, MET_CellOut_Eflow_ety, MET_CellOut_Eflow_sumet)

        # These set up the systematic "SoftTerms_ptHard"
        self.systUtil.setNvtx(nvtxsoftmet)
        #  if(!isData)
        #    self.systUtil.setMETTerm(METUtil::Truth, MET_Truth_NonInt_etx, MET_Truth_NonInt_ety, MET_Truth_NonInt_sumet)

        self.systUtil.setObjectEnergyUncertainties(METUtil::Jets, jesUp, jesDown)
        self.systUtil.setObjectResolutionShift(METUtil::Jets, jerUp, jerDown)

        self.systUtil.setObjectEnergyUncertainties(METUtil::Electrons, eesUp, eesDown)
        self.systUtil.setObjectResolutionShift(METUtil::Electrons, eerUp, eerDown)

        self.systUtil.setObjectEnergyUncertainties(METUtil::Photons, pesUp, pesDown)
        self.systUtil.setObjectResolutionShift(METUtil::Photons, perUp, perDown)

        # Muons are complicated, and MET makes use of track, spectro, and combined quantites,
        # so we need all of their resolutions.
        # comboms reflects that it is the combined muon res affected by shifting ms res up and down.
        # comboid is for shifting the id res up and down
        self.systUtil.setObjectResolutionShift(METUtil::MuonsComboMS, cb_mermsUp, cb_mermsDown)
        self.systUtil.setObjectResolutionShift(METUtil::MuonsComboID, cb_meridUp, cb_meridDown)
        self.systUtil.setObjectResolutionShift(METUtil::SpectroMuons, mermsUp, mermsDown)

        # For now the mes only affects combined
        self.systUtil.setObjectEnergyUncertainties(METUtil::Muons, mesUp, mesDown)

        # Taus are just in as an example
        self.systUtil.setObjectEnergyUncertainties(METUtil::Taus, tesUp, tesDown)

        # Fill histograms

        METObject met_RefFinal = self.systUtil.getMissingET(METUtil::RefFinal)
        h_RefFinal.Fill(met_RefFinal.et()/1000)
        # Jet systematics
        METObject met_RefFinal_JESUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::JESUp)
        h_RefFinal_JESUp.Fill(met_RefFinal_JESUp.et()/1000)
        METObject met_RefFinal_JESDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::JESDown)
        h_RefFinal_JESDown.Fill(met_RefFinal_JESDown.et()/1000)
        METObject met_RefFinal_JERUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::JERUp)
        h_RefFinal_JERUp.Fill(met_RefFinal_JERUp.et()/1000)
        # Electron systematics
        METObject met_RefFinal_EESUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::EESUp)
        h_RefFinal_EESUp.Fill(met_RefFinal_EESUp.et()/1000)
        METObject met_RefFinal_EESDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::EESDown)
        h_RefFinal_EESDown.Fill(met_RefFinal_EESDown.et()/1000)
        METObject met_RefFinal_EERUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::EERUp)
        h_RefFinal_EERUp.Fill(met_RefFinal_EERUp.et()/1000)
        METObject met_RefFinal_EERDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::EERDown)
        h_RefFinal_EERDown.Fill(met_RefFinal_EERDown.et()/1000)
        # Photon systematics
        METObject met_RefFinal_PESUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::PESUp)
        h_RefFinal_PESUp.Fill(met_RefFinal_PESUp.et()/1000)
        METObject met_RefFinal_PESDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::PESDown)
        h_RefFinal_PESDown.Fill(met_RefFinal_PESDown.et()/1000)
        METObject met_RefFinal_PERUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::PERUp)
        h_RefFinal_PERUp.Fill(met_RefFinal_PERUp.et()/1000)
        METObject met_RefFinal_PERDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::PERDown)
        h_RefFinal_PERDown.Fill(met_RefFinal_PERDown.et()/1000)
        # Muon systematics
        METObject met_RefFinal_MESUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::MESUp)
        h_RefFinal_MESUp.Fill(met_RefFinal_MESUp.et()/1000)
        METObject met_RefFinal_MESDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::MESDown)
        h_RefFinal_MESDown.Fill(met_RefFinal_MESDown.et()/1000)
        METObject met_RefFinal_MERIDUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::MERIDUp)
        h_RefFinal_MERIDUp.Fill(met_RefFinal_MERIDUp.et()/1000)
        METObject met_RefFinal_MERIDDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::MERIDDown)
        h_RefFinal_MERIDDown.Fill(met_RefFinal_MERIDDown.et()/1000)
        METObject met_RefFinal_MERMSUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::MERMSUp)
        h_RefFinal_MERMSUp.Fill(met_RefFinal_MERMSUp.et()/1000)
        METObject met_RefFinal_MERMSDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::MERMSDown)
        h_RefFinal_MERMSDown.Fill(met_RefFinal_MERMSDown.et()/1000)

        # ResoSoftTerms uses gRandom for smearing. Set the seed here however you like.
        if(isData) gRandom.SetSeed(UInt_t(RunNumber * EventNumber))
        else gRandom.SetSeed(UInt_t(mc_channel_number * EventNumber))
        # Soft terms systematics
        METObject met_RefFinal_ScaleSoftTermsUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsUp)
        h_RefFinal_ScaleSoftTermsUp.Fill(met_RefFinal_ScaleSoftTermsUp.et()/1000)
        METObject met_RefFinal_ScaleSoftTermsDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsDown)
        h_RefFinal_ScaleSoftTermsDown.Fill(met_RefFinal_ScaleSoftTermsDown.et()/1000)
        METObject met_RefFinal_ResoSoftTermsUp = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUp)
        h_RefFinal_ResoSoftTermsUp.Fill(met_RefFinal_ResoSoftTermsUp.et()/1000)
        METObject met_RefFinal_ResoSoftTermsDown = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsDown)
        h_RefFinal_ResoSoftTermsDown.Fill(met_RefFinal_ResoSoftTermsDown.et()/1000)
        METObject met_RefFinal_ScaleSoftTermsUp_ptHard = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsUp_ptHard)
        # Alternate soft terms systematics
        h_RefFinal_ScaleSoftTermsUp_ptHard.Fill(met_RefFinal_ScaleSoftTermsUp_ptHard.et()/1000)
        METObject met_RefFinal_ScaleSoftTermsDown_ptHard = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsDown_ptHard)
        h_RefFinal_ScaleSoftTermsDown_ptHard.Fill(met_RefFinal_ScaleSoftTermsDown_ptHard.et()/1000)
        METObject met_RefFinal_ResoSoftTermsUp_ptHard = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUp_ptHard)
        h_RefFinal_ResoSoftTermsUp_ptHard.Fill(met_RefFinal_ResoSoftTermsUp_ptHard.et()/1000)
        METObject met_RefFinal_ResoSoftTermsUpDown_ptHard = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUpDown_ptHard)
        h_RefFinal_ResoSoftTermsUpDown_ptHard.Fill(met_RefFinal_ResoSoftTermsUpDown_ptHard.et()/1000)
        METObject met_RefFinal_ResoSoftTermsDownUp_ptHard = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsDownUp_ptHard)
        h_RefFinal_ResoSoftTermsDownUp_ptHard.Fill(met_RefFinal_ResoSoftTermsDownUp_ptHard.et()/1000)
        METObject met_RefFinal_ResoSoftTermsDown_ptHard = self.systUtil.getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsDown_ptHard)
        h_RefFinal_ResoSoftTermsDown_ptHard.Fill(met_RefFinal_ResoSoftTermsDown_ptHard.et()/1000)

        # Print out some information on each event
        if self.verbose:
            print >> self.stream, "Demonstration of smearing and systematics"
            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "All these are the scalar MET term"
            print >> self.stream, "RefEle (smeared): " << self.systUtil.getMissingET(METUtil::RefEle).et()
            #     print >> self.stream, "RefGamma: " << self.systUtil.getMissingET(METUtil::RefGamma).et()
            print >> self.stream, "RefTau: " << self.systUtil.getMissingET(METUtil::RefTau).et()
            print >> self.stream, "RefJet: " << self.systUtil.getMissingET(METUtil::RefJet).et()
            print >> self.stream, "SoftJets: " << self.systUtil.getMissingET(METUtil::SoftJets).et()
            print >> self.stream, "RefMuon: " << self.systUtil.getMissingET(METUtil::RefMuon).et()
            print >> self.stream, "MuonBoy (smeared and scaled): " << self.systUtil.getMissingET(METUtil::MuonTotal).et()
            print >> self.stream, "CellOut_Eflow: " << self.systUtil.getMissingET(METUtil::CellOutEflow).et()
            print >> self.stream, "RefFinal: " << self.systUtil.getMissingET(METUtil::RefFinal).et()
            print >> self.stream, "Truth: " << self.systUtil.getMissingET(METUtil::Truth).et()
            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "HardTerms: " << self.systUtil.getMissingET(METUtil::HardTerms).et()
            print >> self.stream, "HardTerms stands for the sum of Ref* and MuonBoy."
            print >> self.stream, "SoftTerms: " << self.systUtil.getMissingET(METUtil::SoftTerms).et()
            print >> self.stream, "SoftTerms stands for the sum of CellOut(_Eflow) and SoftJets."
            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream,

            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "RefJet JESUp: " << self.systUtil.getMissingET(METUtil::RefJet,METUtil::JESUp).et()
             << ",  JESDown: " << self.systUtil.getMissingET(METUtil::RefJet,METUtil::JESDown).et()
            print >> self.stream, "RefJet JES Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefJet,METUtil::JES).et()
            print >> self.stream, "RefFinal JESUp: " << met_RefFinal_JESUp.et()
             << ", JESDown: " << met_RefFinal_JESDown.et()
            print >> self.stream, "RefFinal JES Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::JES).et()
            print >> self.stream, "RefJet JERUp: " << self.systUtil.getMissingET(METUtil::RefJet,METUtil::JERUp).et()
             << ", JERDown: " << self.systUtil.getMissingET(METUtil::RefJet,METUtil::JERDown).et()
            print >> self.stream, "RefJet JER Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefJet,METUtil::JER).et()
            print >> self.stream, "RefFinal JERUp: " << met_RefFinal_JERUp.et()
             << ", JERDown: " << self.systUtil.getMissingET(METUtil::RefFinal,METUtil::JERDown).et()
            print >> self.stream, "RefFinal JER Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::JER).et()
            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "RefEle EESUp: " << self.systUtil.getMissingET(METUtil::RefEle,METUtil::EESUp).et()
             << ",  EESDown: " << self.systUtil.getMissingET(METUtil::RefEle,METUtil::EESDown).et()
            print >> self.stream, "RefFinal EESUp: " << met_RefFinal_EESUp.et()
             << ",  EESDown: " << met_RefFinal_EESDown.et()
            print >> self.stream, "RefFinal EES Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::EES).et()
            print >> self.stream, "RefEle EERUp: " << self.systUtil.getMissingET(METUtil::RefEle,METUtil::EERUp).et()
             << ",  EERDown: " << self.systUtil.getMissingET(METUtil::RefEle,METUtil::EERDown).et()
            print >> self.stream, "RefFinal EERUp: " << met_RefFinal_EERUp.et()
             << ",  EERDown: " << met_RefFinal_EERDown.et()
            print >> self.stream, "RefFinal EER Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::EER).et()
            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "MuonBoy MESUp: " << self.systUtil.getMissingET(METUtil::MuonTotal,METUtil::MESUp).et()
             << ",  MESDown: " << self.systUtil.getMissingET(METUtil::MuonTotal,METUtil::MESDown).et()
            print >> self.stream, "RefFinal MESUp: " << met_RefFinal_MESUp.et()
             << ",  MESDown: " << met_RefFinal_MESDown.et()
            print >> self.stream, "RefFinal MES Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::MES).et()
            print >> self.stream, "MuonBoy MERIDUp: " << self.systUtil.getMissingET(METUtil::MuonTotal,METUtil::MERIDUp).et()
             << ",  MERIDDown: " << self.systUtil.getMissingET(METUtil::MuonTotal,METUtil::MERIDDown).et()
            print >> self.stream, "RefFinal MERIDUp: " << met_RefFinal_MERIDUp.et()
             << ",  MERIDDown: " << met_RefFinal_MERIDDown.et()
            print >> self.stream, "RefFinal MERID Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::MERID).et()
            print >> self.stream, "MuonBoy MERMSUp: " << self.systUtil.getMissingET(METUtil::MuonTotal,METUtil::MERMSUp).et()
             << ",  MERMSDown: " << self.systUtil.getMissingET(METUtil::MuonTotal,METUtil::MERMSDown).et()
            print >> self.stream, "RefFinal MERMSUp: " << met_RefFinal_MERMSUp.et()
             << ",  MERMSDown: " << met_RefFinal_MERMSDown.et()
            print >> self.stream, "RefFinal MERMS Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::MERMS).et()
            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "RefTau TESUp: " << self.systUtil.getMissingET(METUtil::RefTau,METUtil::TESUp).et()
             << ",  TESDown: " << self.systUtil.getMissingET(METUtil::RefTau,METUtil::TESDown).et()
            print >> self.stream, "RefFinal TESUp: " << self.systUtil.getMissingET(METUtil::RefFinal,METUtil::TESUp).et()
             << ",  TESDown: " << self.systUtil.getMissingET(METUtil::RefFinal,METUtil::TESDown).et()
            print >> self.stream, "RefFinal TES Diff (up - Down)/none : " << self.systUtil.getMissingETDiff(METUtil::RefFinal,METUtil::TES).et()
            print >> self.stream, "+++++++++++++++++++++++++++++"


            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "Now follow the soft terms. Apart from AllClusters, these already include pileup contributions."
            print >> self.stream, endl
            print >> self.stream, "AllClusters is the PLHC systematic on CellOut(_Eflow) and SoftJets."
            print >> self.stream, "RefFinal AllClusters Up: " << self.systUtil.getMissingET(METUtil::RefFinal,METUtil::AllClustersUp).et()
             << ", RefFinal AllClusters Down: " << self.systUtil.getMissingET(METUtil::RefFinal,METUtil::AllClustersDown).et()

            print >> self.stream, endl
            print >> self.stream, "These are the April 2012 systematics. For information, please see:"
            print >> self.stream, "https:#indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=161247"
            print >> self.stream, endl

            print >> self.stream, "ScaleSoftTerms is the systematic on the scale CellOut(_Eflow) and SoftJets."

                # ResoSoftTerms uses gRandom for smearing. Set the seed here however you like.
            if(isData) gRandom.SetSeed(UInt_t(RunNumber * EventNumber))
            else gRandom.SetSeed(UInt_t(mc_channel_number * EventNumber))

            print >> self.stream, "RefFinal ScaleSoftTerms Up: " << met_RefFinal_ScaleSoftTermsUp.et()
             << ", RefFinal ScaleSoftTerms Down: " << met_RefFinal_ScaleSoftTermsDown.et()
            print >> self.stream, "ResoSoftTerms is the systematic on the scale CellOut(_Eflow) and SoftJets."
            print >> self.stream, "RefFinal ResoSoftTerms Up: " << met_RefFinal_ResoSoftTermsUp.et()
             << ", RefFinal ResoSoftTerms Down: " << met_RefFinal_ResoSoftTermsDown.et()

            print >> self.stream, "ScaleSoftTerms_ptHard is the systematic on the scale CellOut(_Eflow) and SoftJets."
            print >> self.stream, "RefFinal ScaleSoftTerms_ptHard Up: " << met_RefFinal_ScaleSoftTermsUp_ptHard.et()
             << ", RefFinal ScaleSoftTerms_ptHard Down: " << met_RefFinal_ScaleSoftTermsDown_ptHard.et()
            print >> self.stream, "ResoSoftTerms is the systematic on the scale CellOut(_Eflow) and SoftJets."
            print >> self.stream, "This is parameterised in terms of longitudinal and transverse components, which can be varied coherently or anti-coherently."
            print >> self.stream, "RefFinal ResoSoftTerms_ptHard Up: " << met_RefFinal_ResoSoftTermsUp_ptHard.et()
             << ", RefFinal ResoSoftTerms_ptHard Down: " << met_RefFinal_ResoSoftTermsDown_ptHard.et()
            print >> self.stream, "RefFinal ResoSoftTerms_ptHard UpDown: " << met_RefFinal_ResoSoftTermsUpDown_ptHard.et()
             << ", RefFinal ResoSoftTerms_ptHard DownUp: " << met_RefFinal_ResoSoftTermsDownUp_ptHard.et()

            print >> self.stream, "+++++++++++++++++++++++++++++"
            print >> self.stream, "Combined errors, giving an uncertainty on MET"
            METObject smearedMET = self.systUtil.getMissingET(METUtil::RefFinal)
            print >> self.stream, "RefFinal MET = " << smearedMET.et() << " +- " << self.systUtil.absDeltaMissingET(METUtil::RefFinal).et()
             << " (" << 100*self.systUtil.deltaMissingET(METUtil::RefFinal).et() << "%)"
            print >> self.stream, "+++++++++++++++++++++++++++++"

    def Terminate(self):

        self.outfile.Write()
        self.outfile.Close()

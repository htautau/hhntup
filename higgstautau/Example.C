#define EventReader_cxx
#include "Example.h"
#include "TMath.h"
#include "TH1D.h"
#include "TFile.h"
#include "../MissingETUtility/METUtility.h"
#include "JetUncertainties/MultijetJESUncertaintyProvider.h"
#include "JetResolution/JERProvider.h"
#include "egammaAnalysisUtils/EnergyRescaler.h"
#include "MuonMomentumCorrections/SmearingClass.h"
#include "../MissingETUtility/TESUncertaintyProvider.h"

void Example::Begin(TTree */*tree*/)
{
    // Initialise your METUtility object
    m_testUtil = new METUtility;

    // *** Demonstration of the configuration here  ***
    // *** All values that are set are the defaults ***

    // Turn on (off) the relevant MET terms
    // Standard MET_RefFinal has:
    // RefEle, RefGamma, RefTau, RefJet, SoftJets, RefMuon, MuonTotal, (CellOut), CellOut_Eflow
    m_testUtil->defineMissingET(true, true, true, true, true, true, true, false, true);

    // SUSY group, MET_Simplified20 has
    // RefEle, (RefGamma), (RefTau), RefJet, (SoftJets), (RefMuon), MuonTotal, CellOut, (CellOut_Eflow)
    // m_testUtil->defineMissingET(true, false, true, true, false, false, true, true, false);

    // The threshold below which jets enter the SoftJets term (JES is not applied)
    m_testUtil->setSoftJetCut(20e3);

    // Whether to use MUID muons (otherwise STACO).
    m_testUtil->setIsMuid(false);

    // Whether METUtility should scream at you over every little thing
    m_testUtil->setVerbosity(false);

    m_systUtil = new METUtility;

    // Some other options are available, but not recommended/deprecated

    // *** Set up the uncertainty tools ***
    // Tag assumed: JetUncertainties-00-05-09-02
    m_jesTool = new MultijetJESUncertaintyProvider("MultijetJES_Preliminary.config",
            "InsituJES2011_AllNuisanceParameters.config",
            "AntiKt4LCTopoJets","MC11c");

    // Tag assumed: JetResolution-01-00-00
    m_jerTool = new JERProvider("AntiKt4LCTopoJES", "Truth", "../../Jet/JetResolution/share/JERProviderPlots.root");
    m_jerTool->init();

    // Tag assumed: egammaAnalysisUtils-00-02-76
    m_egammaTool = new eg2011::EnergyRescaler;
    m_egammaTool->useDefaultCalibConstants("2011");

    // Tag assumed: MuonMomentumCorrections-00-05-03
    m_muonTool = new MuonSmear::SmearingClass("Data11","staco","q_pT","Rel17");
    m_muonTool->UseScale(1);
    m_muonTool->UseImprovedCombine();

    // No tag yet, testing code
    m_tesTool = new TESUncertaintyProvider;

    // Prepare an output file
    m_outfile = new TFile("metutil_example.root","RECREATE");
    // Set up a bunch of histograms
    h_RefFinal = new TH1D("h_RefFinal",
            "RefFinal",50,0,1000);
    // Jet systematics
    h_RefFinal_JESUp = new TH1D("h_RefFinal_JESUp",
            "RefFinal JESUp",50,0,1000);
    h_RefFinal_JESDown = new TH1D("h_RefFinal_JESDown",
            "RefFinal JESDown",50,0,1000);
    h_RefFinal_JERUp = new TH1D("h_RefFinal_JERUp",
            "RefFinal JERUp",50,0,1000);
    // Electron systematics
    h_RefFinal_EESUp = new TH1D("h_RefFinal_EESUp",
            "RefFinal EESUp",50,0,1000);
    h_RefFinal_EESDown = new TH1D("h_RefFinal_EESDown",
            "RefFinal EESDown",50,0,1000);
    h_RefFinal_EERUp = new TH1D("h_RefFinal_EERUp",
            "RefFinal EERUp",50,0,1000);
    h_RefFinal_EERDown = new TH1D("h_RefFinal_EERDown",
            "RefFinal EERDown",50,0,1000);
    // Photon systematics
    h_RefFinal_PESUp = new TH1D("h_RefFinal_PESUp",
            "RefFinal PESUp",50,0,1000);
    h_RefFinal_PESDown = new TH1D("h_RefFinal_PESDown",
            "RefFinal PESDown",50,0,1000);
    h_RefFinal_PERUp = new TH1D("h_RefFinal_PERUp",
            "RefFinal PERUp",50,0,1000);
    h_RefFinal_PERDown = new TH1D("h_RefFinal_PERDown",
            "RefFinal PERDown",50,0,1000);
    // Muon systematics
    h_RefFinal_MESUp = new TH1D("h_RefFinal_MESUp",
            "RefFinal MESUp",50,0,1000);
    h_RefFinal_MESDown = new TH1D("h_RefFinal_MESDown",
            "RefFinal MESDown",50,0,1000);
    h_RefFinal_MERIDUp = new TH1D("h_RefFinal_MERIDUp",
            "RefFinal MERIDUp",50,0,1000);
    h_RefFinal_MERIDDown = new TH1D("h_RefFinal_MERIDDown",
            "RefFinal MERIDDown",50,0,1000);
    h_RefFinal_MERMSUp = new TH1D("h_RefFinal_MERMSUp",
            "RefFinal MERMSUp",50,0,1000);
    h_RefFinal_MERMSDown = new TH1D("h_RefFinal_MERMSDown",
            "RefFinal MERMSDown",50,0,1000);
    // Soft terms systematics
    h_RefFinal_ScaleSoftTermsUp = new TH1D("h_RefFinal_ScaleSoftTermsUp",
            "RefFinal ScaleSoftTermsUp",50,0,1000);
    h_RefFinal_ScaleSoftTermsDown = new TH1D("h_RefFinal_ScaleSoftTermsDown",
            "RefFinal ScaleSoftTermsDown",50,0,1000);
    h_RefFinal_ResoSoftTermsUp = new TH1D("h_RefFinal_ResoSoftTermsUp",
            "RefFinal ResoSoftTermsUp",50,0,1000);
    h_RefFinal_ResoSoftTermsDown = new TH1D("h_RefFinal_ResoSoftTermsDown",
            "RefFinal ResoSoftTermsDown",50,0,1000);
    // Alternate soft terms systematics
    h_RefFinal_ScaleSoftTermsUp_ptHard = new TH1D("h_RefFinal_ScaleSoftTermsUp_ptHard",
            "RefFinal ScaleSoftTermsUp_ptHard",50,0,1000);
    h_RefFinal_ScaleSoftTermsDown_ptHard = new TH1D("h_RefFinal_ScaleSoftTermsDown_ptHard",
            "RefFinal ScaleSoftTermsDown_ptHard",50,0,1000);
    h_RefFinal_ResoSoftTermsUp_ptHard = new TH1D("h_RefFinal_ResoSoftTermsUp_ptHard",
            "RefFinal ResoSoftTermsUp_ptHard",50,0,1000);
    h_RefFinal_ResoSoftTermsUpDown_ptHard = new TH1D("h_RefFinal_ResoSoftTermsUpDown_ptHard",
            "RefFinal ResoSoftTermsUpDown_ptHard",50,0,1000);
    h_RefFinal_ResoSoftTermsDownUp_ptHard = new TH1D("h_RefFinal_ResoSoftTermsDownUp_ptHard",
            "RefFinal ResoSoftTermsDownUp_ptHard",50,0,1000);
    h_RefFinal_ResoSoftTermsDown_ptHard = new TH1D("h_RefFinal_ResoSoftTermsDown_ptHard",
            "RefFinal ResoSoftTermsDown_ptHard",50,0,1000);
}

Bool_t Example::Process(Long64_t entry)
{
    bool verbose = true;
    bool isData = false;

    // This just makes sure aliases are set to some generic objects
    EventReader::Process(entry);

    // See this method for an example of setting up the METUtility
    // such that it can rebuild the D3PD nominal MET values.
    Example::CheckMetRebuilding(verbose);

    // See this method for an example of how to pass object systematics
    // to METUtility such that they can be propagated to the MET.
    Example::DemoMetSystematics(verbose,isData);

    return kTRUE;
}

void Example::CheckMetRebuilding(bool verbose)
{

    //////////////////////////////////////////////////////////////////////////////
    // Demonstrates how to set up the METUtility with object momenta such that
    // MET_RefFinal can be rebuilt matching the values in D3PD.

    // Start with a clean METUtility
    m_testUtil->reset();

    // For these, just the kinematics need to be set
    m_testUtil->setElectronParameters(el_pt, el_eta, el_phi,
            el_MET_wet, el_MET_wpx, el_MET_wpy,
            el_MET_statusWord);

    m_testUtil->setPhotonParameters(ph_pt, ph_eta, ph_phi,
            ph_MET_wet, ph_MET_wpx, ph_MET_wpy,
            ph_MET_statusWord);

    //  Tau rebuilding is unsafe. The tau finding is frequently rerun in D3PD.
    //  m_testUtil->setTauParameters(tau_pt, tau_eta, tau_phi,
    //			     tau_MET_wet, tau_MET_wpx, tau_MET_wpy,
    //			     tau_MET_statusWord);

    //  Cluster rebuilding is just for specialist studies
    //  m_testUtil->setClusterParameters(cl_pt, cl_eta, cl_phi,
    //				 cl_MET_wet, cl_MET_wpx, cl_MET_wpy,
    //				 cl_MET_statusWord);

    //  Track rebuilding is just for specialist studies
    //  m_testUtil->setTrackParameters(trk_pt, trk_eta, trk_phi,
    //			       trk_MET_wet, trk_MET_wpx, trk_MET_wpy,
    //			       trk_MET_statusWord);

    // The SoftJets term is now to be taken from D3PD, and no "hard jets" are allowed
    // to enter it. Recalibration or smearing could cause jets that were above the
    // 20 GeV threshold to drop below it, so we supply the original pT's to prevent
    // them from being moved out of RefJet.
    m_testUtil->setJetParameters(jet_pt, jet_eta, jet_phi, jet_E,
            jet_MET_wet, jet_MET_wpx, jet_MET_wpy,
            jet_MET_statusWord);
    m_testUtil->setOriJetParameters(jet_pt);

    // Muons may be ID, combined, or standalone. For the latter especially,
    // we need to set the MS four-momenta because they are treated differently.
    m_testUtil->setMuonParameters(mu_pt, mu_eta, mu_phi,
            mu_MET_wet, mu_MET_wpx, mu_MET_wpy,
            mu_MET_statusWord);
    m_testUtil->setExtraMuonParameters(mu_ms_qoverp, mu_ms_theta, mu_ms_phi, mu_charge);
    // An alternative version of this method is demonstrated below, and takes pT, eta, phi instead.
    // This is more convenient when one needs to smear the pT, for example.

    // When the terms are not rebuilt from the objects, due to incomplete info,
    // then one needs to set the term directly from the stored value in the D3PD.
    // This might also be done if you aren't interested in the possible variations
    // of that term. E.g. if you only care about photons, no need to rebuild muons.

    // m_testUtil->setMETTerm(METUtil::RefJet, MET_RefJet_etx, MET_RefJet_ety, MET_RefJet_sumet);  
    m_testUtil->setMETTerm(METUtil::SoftJets, MET_SoftJets_etx, MET_SoftJets_ety, MET_SoftJets_sumet);  
    // m_testUtil->setMETTerm(METUtil::RefEle, MET_RefEle_etx, MET_RefEle_ety, MET_RefEle_sumet);  
    // m_testUtil->setMETTerm(METUtil::RefGamma, MET_RefGamma_etx, MET_RefGamma_ety, MET_RefGamma_sumet);

    // *** Note the difference in naming -- there is normally no MET_MuonTotal term.
    //     It's usually Muon_Total, Muon_Total_Muid, something like that.
    //     MET_RefFinal in particular uses MET_MuonBoy.
    // m_testUtil->setMETTerm(METUtil::MuonTotal, MET_MuonBoy_etx, MET_MuonBoy_ety, MET_MuonBoy_sumet);

    // *** Note that RefMuon is not rebuilt from muons -- it is a calorimeter term.
    m_testUtil->setMETTerm(METUtil::RefMuon, MET_RefMuon_etx, MET_RefMuon_ety, MET_RefMuon_sumet);
    m_testUtil->setMETTerm(METUtil::RefTau, MET_RefTau_etx, MET_RefTau_ety, MET_RefTau_sumet);
    // m_testUtil->setMETTerm(METUtil::CellOut, MET_CellOut_etx, MET_CellOut_ety, MET_CellOut_sumet);
    m_testUtil->setMETTerm(METUtil::CellOutEflow, MET_CellOut_Eflow_etx, MET_CellOut_Eflow_ety, MET_CellOut_Eflow_sumet);

    // This is the simple check, where you compare the terms manually against what's in the D3PD.
    // Note: every call to getMissingET recomputes the terms, so if you need to get more than one
    // value, e.g. etx, ety, et, sumET, it's more efficient to get the METObject.
    // Usually, comparing etx and/or ety is more informative, because et could be right if
    // etx and ety were flipped, for example. They also add linearly, which et doesn't.

    METObject RefEle_util = m_testUtil->getMissingET(METUtil::RefEle);
    METObject RefGamma_util = m_testUtil->getMissingET(METUtil::RefGamma);
    METObject RefTau_util = m_testUtil->getMissingET(METUtil::RefTau);
    METObject RefMuon_util = m_testUtil->getMissingET(METUtil::RefMuon);
    METObject RefJet_util = m_testUtil->getMissingET(METUtil::RefJet);
    METObject SoftJets_util = m_testUtil->getMissingET(METUtil::SoftJets);
    //      METObject refCellOut_util = m_testUtil->getMissingET(METUtil::CellOut);
    METObject CellOutEflow_util = m_testUtil->getMissingET(METUtil::CellOutEflow);
    METObject MuonTotal_util = m_testUtil->getMissingET(METUtil::MuonTotal);
    METObject RefFinal_util = m_testUtil->getMissingET(METUtil::RefFinal);

    if(verbose)  {
        cout << "** Manual consistency check **" << endl << endl;
        cout << "Term:    Original   vs    Tool output" << endl;
        cout << "RefEle etx: "    << MET_RefEle_etx   << " vs " << RefEle_util.etx()    << endl;
        cout << "RefGamma etx: "  << MET_RefGamma_etx << " vs " << RefGamma_util.etx()  << endl;
        cout << "RefTau etx: "    << MET_RefTau_etx   << " vs " << RefTau_util.etx()    << endl;
        cout << "RefMuon etx: "   << MET_RefMuon_etx  << " vs " << RefMuon_util.etx()   << endl;
        cout << "RefJet etx: "    << MET_RefJet_etx   << " vs " << RefJet_util.etx()   << endl;
        cout << "SoftJets etx: "   << MET_SoftJets_etx << " vs " << SoftJets_util.etx()   << endl;
        cout << "MuonBoy etx: "   << MET_MuonBoy_etx  << " vs " << MuonTotal_util.etx() << endl;
        cout << "CellOut_Eflow etx: " << MET_CellOut_Eflow_etx << " vs " << CellOutEflow_util.etx() << endl;
        cout << "RefFinal etx: "  << MET_RefFinal_etx << " vs " << RefFinal_util.etx()  << endl << endl;
    }

    // If you don't want to test manually, there's a built-in consistency check.
    // To test just one term, fill a METObject with the appropriate values,
    // then feed it to the checkConsistency() method.
    // The difference can be retrieved via a reference argument.

    METObject refFinal_test(MET_RefFinal_etx,
            MET_RefFinal_ety,
            MET_RefFinal_sumet);
    bool check_refFinal = m_testUtil->checkConsistency(refFinal_test,METUtil::RefFinal);
    if(check_refFinal) cout << "RefFinal checks out!" << endl;
    else cout << "RefFinal doesn't check out!" << endl;

    // By filling a vector of terms, you can test all of them in one function call.
    // The sum (and sumET) will be tested as well. (Can't get the difference this way).
    METObject refEle_test(MET_RefEle_etx,
            MET_RefEle_ety,
            MET_RefEle_sumet);
    METObject refGamma_test(MET_RefGamma_etx,
            MET_RefGamma_ety,
            MET_RefGamma_sumet);
    METObject refJet_test(MET_RefJet_etx,
            MET_RefJet_ety,
            MET_RefJet_sumet);
    METObject muonBoy_test(MET_MuonBoy_etx,
            MET_MuonBoy_ety,
            MET_MuonBoy_sumet);
    vector<pair<int,METObject> > testvector;
    testvector.push_back(pair<int,METObject>(METUtil::RefEle,refEle_test));
    testvector.push_back(pair<int,METObject>(METUtil::RefGamma,refGamma_test));
    testvector.push_back(pair<int,METObject>(METUtil::RefJet,refJet_test));
    testvector.push_back(pair<int,METObject>(METUtil::MuonTotal,muonBoy_test));

    bool check = m_testUtil->checkConsistency(testvector);
    if(check) cout << "MET checks out!" << endl;
    else cout << "MET doesn't check out!" << endl;

    // In addition to the etx, ety, sumet retrieval, you can also get the MET significance.
    // By default, METObject::sig() returns etx() / ( 0.5*sqrt(sumet()) )
    // 
    // There is also the possibility of returning a more sophisticated estimator for
    // the significance, activated by calling METUtility::doSignificance() in setup.
    // This is still under development, and requires all object resolutions to be set.

}

void Example::DemoMetSystematics(bool verbose,bool isData)
{

    //////////////////////////////////////////////////////////////////////////////
    // Demonstrates how to set up the METUtility with object momenta
    // such that MET_RefFinal can be rebuilt matching the values in D3PD.
    //
    // *** *** *** *** *** DISCLAIMER *** *** *** *** ***
    //
    // These examples of uncertainty-setting are meant to
    // demonstrate how the uncertainties should be passed
    // to METUtility.  Recommendations on just which ones
    // you are meant to be using come from the CP groups.

    // Check for a good primary vertex
    // This is needed for jet and soft term systematics
    bool goodPV = false;
    int nvtxsoftmet = 0;
    int nvtxjets = 0;
    if(vxp_n>0) {
        // Most D3PDs contain the vx_type branch, but some don't.
        // Those which don't are most likely skimmed, and require at least 1 primary vertex for all events.
        // If your D3PD is skimmed in this way, then the goodPV (nTracks and z) check should be applied
        // to the first vertex in the collection.
        // Otherwise, you should ensure that the vx_type branch is available.
        for(int i=0; i<vxp_n; i++) {
            if(vxp_type->at(i) == 1 && vxp_nTracks->at(i)>2 && fabs(vxp_z->at(i))<200.) goodPV = true;
        }
        if(goodPV) {
            for(int i=0; i<vxp_n; i++) {
                if(vxp_nTracks->at(i)>2) nvtxsoftmet++;
                if(vxp_nTracks->at(i)>1) nvtxjets++;
            }
        }
    }  
    // First, we get the jet energy scale uncertainties and
    // systematic variation in the jet resolutions
    vector<float> jesUp;
    vector<float> jesDown;
    vector<float> jerUp;
    vector<float> jerDown;

    // Note on use of ROOT random number generators:
    // TRandom and TRandom2 have many documented deficiencies.
    // TRandom3 is generally considered safely usable.
    // Also note that ROOT's gRandom calls TRandom3.
    TRandom3 *jetRandom = new TRandom3;

    for(int iJet = 0; iJet < jet_n; ++iJet){
        float jesShiftUp = 0.0;
        float jesShiftDown = 0.0;
        float jerShift = 1.0;
        // Safest to assume nothing about the uncertainties on soft jets.
        // These will go into SoftJets anyhow, and so the JES systematics
        // aren't used.

        if(jet_pt->at(iJet) > 20e3
                && jet_pt->at(iJet) < 7000e3
                && fabs(jet_eta->at(iJet)) < 4.5){

            // delta R cut needed to apply close-by jets uncertainty
            float drmin=9999;
            double pi = TMath::Pi();
            if( jet_pt->at(iJet)>20000) {
                for (int ii = 0; ii < jet_n; ii++ ) {
                    if(jet_emscale_pt->at(ii)>7000) {
                        if(iJet!=ii) {
                            double deta = jet_eta->at(iJet) - jet_eta->at(ii);
                            double dphi = fabs(fmod((jet_phi->at(iJet)
                                            - jet_phi->at(ii))+3*pi,2*pi)-pi);
                            double dr = sqrt(deta*deta+dphi*dphi);
                            if(dr<drmin) drmin=dr;
                        }
                    }
                }
            }

            // The bool is the "isPos" argument
            jesShiftUp = m_jesTool->getRelUncert(jet_pt->at(iJet),
                    jet_eta->at(iJet),drmin,
                    true, nvtxjets, averageIntPerXing);
            jesShiftDown = -1*m_jesTool->getRelUncert(jet_pt->at(iJet),
                    jet_eta->at(iJet),drmin,
                    false, nvtxjets, averageIntPerXing);
        }
        jesUp.push_back(jesShiftUp);
        jesDown.push_back(jesShiftDown);

        // Allowable range is > 10 GeV, but anything below 20 enters SoftJets
        if(jet_pt->at(iJet) > 20e3 && jet_pt->at(iJet) < 5000e3){
            double pt = jet_pt->at(iJet);
            double eta = jet_eta->at(iJet);
            if(fabs(eta)>4.5) eta = eta>0 ? 4.49 : -4.49;

            double S = m_jerTool->getRelResolutionMC(pt/1e3,eta);
            double U = m_jerTool->getResolutionUncert(pt/1e3,eta);
            double smearingFactorSyst = sqrt(pow(S+U,2)-pow(S,2));

            // You can set the seed however you like, but if reproducibility
            // is required, setting it to something like object phi ensures
            // a good mix of randomness and reproducibility.
            jetRandom->SetSeed(int(1.e5*jet_phi->at(iJet)));
            jerShift = jetRandom->Gaus(0, smearingFactorSyst);
        }

        jerUp.push_back(jerShift);
        jerDown.push_back(-1*jerShift); // Usually not used, see below.

        //////////////////////////////////////////////////////////////////////
        // Note: The JERDown shift is essentially meaningless.
        // If one is smearing central values, then there is an alternate
        // definition, i.e. from r16:
        //
        // S = m_jerTool.getRelResolutionData(pt/1e3,eta);
        // SMC = m_jerTool.getRelResolutionMC(pt/1e3,eta);
        // U = m_jerTool.getResolutionUncert(pt/1e3,eta);
        // smearingFactorMC = sqrt( S*S - SMC*SMC );
        // smearingFactorSystUp = sqrt( (S+U)*(S+U) - SMC*SMC );
        // smearingFactorSystDown = (S-U > SMC) ? sqrt( (S+U)*(S+U) - SMC*SMC ) : 0;
        // 
        // float jerShift = jetRandom->Gaus(1,smearingFactorMC);
        // float jerShiftUp = jetRandom->Gaus(1,smearingFactorSystUp)/jerShift;
        // float jerShiftDown = jetRandom->Gaus(1,smearingFactorSystDown)/jerShift;
        //
        // jet_smeared_pt = pt*jerShift;
        // jerUp.push_back(jerShiftUp-1);
        // jerDown.push_back(jerShiftDown-1);
        //
        // This means that we smear the MC jets to match the resolution in data
        // for central values, or the resolution +/- uncertainty.
        // The standard practice is only to use res + uncertainty.
        //
        //////////////////////////////////////////////////////////////////////

    }//end of jet loop

    delete jetRandom;

    // Here we get the electron energy scale and resolution systematics
    vector<float> eesUp;
    vector<float> eesDown;
    vector<float> eerUp;
    vector<float> eerDown;
    vector<float> *el_smeared_pt = new vector<float>;

    for (unsigned int iEl = 0; iEl < el_pt->size(); ++iEl) {

        m_egammaTool->SetRandomSeed(int(1e5*fabs(el_phi->at(iEl))));

        // Smear to match the data resolution, or by systematic variations
        float smear = m_egammaTool->getSmearingCorrectionMeV(el_cl_eta->at(iEl), el_E->at(iEl), 0, true);
        float smearUp = m_egammaTool->getSmearingCorrectionMeV(el_cl_eta->at(iEl), el_E->at(iEl), 2, true);
        float smearDown = m_egammaTool->getSmearingCorrectionMeV(el_cl_eta->at(iEl), el_E->at(iEl), 1, true);

        el_smeared_pt->push_back(smear*el_pt->at(iEl));
        eerUp.push_back((smearUp - smear)/smear);
        eerDown.push_back((smearDown - smear)/smear);

        // Correct the measured energies in data, and scale by systematic variations
        float correction = 1.;
        if(isData)
            correction = m_egammaTool->applyEnergyCorrectionMeV(el_cl_eta->at(iEl),el_cl_phi->at(iEl),
                    el_E->at(iEl),el_cl_pt->at(iEl),0,"ELECTRON") / el_E->at(iEl);
        el_smeared_pt->at(iEl)*= correction;
        double energyUp = m_egammaTool->applyEnergyCorrectionMeV(el_cl_eta->at(iEl),el_cl_phi->at(iEl),
                el_E->at(iEl),el_cl_pt->at(iEl),2,"ELECTRON") / (correction*el_E->at(iEl)) - 1;
        double energyDown = m_egammaTool->applyEnergyCorrectionMeV(el_cl_eta->at(iEl),el_cl_phi->at(iEl),
                el_E->at(iEl),el_cl_pt->at(iEl),1,"ELECTRON") / (correction*el_E->at(iEl)) - 1;

        eesUp.push_back(energyUp);
        eesDown.push_back(energyDown);
    } //end of electron loop


    // Now we get the same for photons
    vector<float> pesUp;
    vector<float> pesDown;
    vector<float> perUp;
    vector<float> perDown;
    vector<float> *ph_smeared_pt = new vector<float>;

    for (unsigned int iPh = 0; iPh < ph_pt->size(); ++iPh) {

        m_egammaTool->SetRandomSeed(int(1.e+5*fabs(ph_phi->at(iPh))));

        // Smear to match the data resolution, or by systematic variations
        float smear = m_egammaTool->getSmearingCorrectionMeV(ph_cl_eta->at(iPh), ph_E->at(iPh), 0, true);
        float smearUp = m_egammaTool->getSmearingCorrectionMeV(ph_cl_eta->at(iPh), ph_E->at(iPh), 2, true);
        float smearDown = m_egammaTool->getSmearingCorrectionMeV(ph_cl_eta->at(iPh), ph_E->at(iPh), 1, true);

        ph_smeared_pt->push_back(smear*ph_pt->at(iPh));
        perUp.push_back((smearUp - smear)/smear);
        perDown.push_back((smearDown - smear)/smear);

        // Correct the measured energies in data, and scale by systematic variations
        // Conversions are treated differently.
        float correction = 1.;
        string photontype = ph_isConv->at(iPh) ? "CONVERTED_PHOTON" : "UNCONVERTED_PHOTON";
        if(isData)
            correction = m_egammaTool->applyEnergyCorrectionMeV(ph_cl_eta->at(iPh),ph_cl_phi->at(iPh),ph_E->at(iPh),ph_cl_pt->at(iPh),0,photontype) / ph_E->at(iPh);
        ph_smeared_pt->at(iPh)*= correction;
        double energyUp = m_egammaTool->applyEnergyCorrectionMeV(ph_cl_eta->at(iPh),ph_cl_phi->at(iPh),
                ph_E->at(iPh),ph_cl_pt->at(iPh),2,photontype) / (correction*ph_E->at(iPh)) - 1;
        double energyDown = m_egammaTool->applyEnergyCorrectionMeV(ph_cl_eta->at(iPh),ph_cl_phi->at(iPh),
                ph_E->at(iPh),ph_cl_pt->at(iPh),1,photontype) / (correction*ph_E->at(iPh)) - 1;

        pesUp.push_back(energyUp);
        pesDown.push_back(energyDown);
    }//end of photon loop

    // And now the same for muons. We need resolution shifts for ID and MS,
    // and different treatment for the MS four-vector (for standalone muons).
    vector<float> *mu_smeared_pt = new vector<float>;
    vector<float> *mu_smeared_ms_pt = new vector<float>;

    vector<float> cb_meridUp;
    vector<float> cb_meridDown;
    vector<float> cb_mermsUp;
    vector<float> cb_mermsDown;
    vector<float> mermsUp;
    vector<float> mermsDown;

    vector<float> mesUp;
    vector<float> mesDown;

    for(unsigned int iMu = 0; iMu < mu_pt->size(); ++iMu){

        double ptcb = mu_pt->at(iMu);
        double ptid = (mu_id_qoverp_exPV->at(iMu) != 0.) ? fabs(sin(mu_id_theta_exPV->at(iMu))/mu_id_qoverp_exPV->at(iMu)) : 0.;
        double ptms = (mu_ms_qoverp->at(iMu) != 0.) ? fabs(sin(mu_ms_theta->at(iMu))/mu_ms_qoverp->at(iMu)) : 0.;
        m_muonTool->SetSeed(int(1.e+5*fabs(mu_phi->at(iMu))));
        double etaMu = mu_eta->at(iMu);
        double charge = mu_charge->at(iMu);
        m_muonTool->Event(ptms, ptid, ptcb, etaMu, charge);

        Float_t smearedCombinedPt = m_muonTool->pTCB();
        if(!mu_isCombinedMuon->at(iMu)) smearedCombinedPt = m_muonTool->pTMS() + m_muonTool->pTID();

        Float_t smearedMSPt = m_muonTool->pTMS();

        mu_smeared_ms_pt->push_back(smearedMSPt);
        mu_smeared_pt->push_back(smearedCombinedPt);

        double ptMS_smeared, ptID_smeared, ptCB_smeared;
        float smearedpTMS, smearedpTID, smearedpTCB;
        smearedpTMS = 0.1; smearedpTID = 0.1; smearedpTCB = 0.1;
        m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSLOW");
        smearedpTMS = ptMS_smeared/smearedMSPt - 1.0;
        smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
        mermsDown.push_back(smearedpTMS);
        cb_mermsDown.push_back(smearedpTCB);
        m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "MSUP");
        smearedpTMS = ptMS_smeared/smearedMSPt - 1.0;
        smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
        mermsUp.push_back(smearedpTMS);
        cb_mermsUp.push_back(smearedpTCB);
        m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDUP");
        smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
        cb_meridUp.push_back(smearedpTCB);
        m_muonTool->PTVar(ptMS_smeared, ptID_smeared, ptCB_smeared, "IDLOW");
        smearedpTCB = ptCB_smeared/smearedCombinedPt - 1.0;
        cb_meridDown.push_back(smearedpTCB);

        int detRegion = m_muonTool->DetRegion();
        if(detRegion==-1) detRegion = 3;
        double scalesyst = m_muonTool->getScaleSyst_CB().at(detRegion);    
        mesUp.push_back(scalesyst);
        mesDown.push_back(-scalesyst);

    }//end of muon loop

    vector<float> tesUp;
    vector<float> tesDown;

    // And for taus (this is test code, do not use without tau group approval)
    for(int iTau=0; iTau<tau_n; iTau++) {
        double pt = tau_pt->at(iTau);
        double eta = tau_eta->at(iTau);
        int nProng = tau_nProng->at(iTau);
        double uncert = m_tesTool->GetTESUncertainty(pt/1e3, eta, nProng);

        if(uncert < 0) uncert = 0;
        tesUp.push_back(uncert);
        tesDown.push_back(-1*uncert);
    }

    // This demonstration is for doing smearing and systematics
    m_systUtil->reset();
    m_systUtil->setJetParameters(jet_pt, jet_eta, jet_phi, jet_E,
            jet_MET_wet, jet_MET_wpx, jet_MET_wpy,
            jet_MET_statusWord);
    m_systUtil->setOriJetParameters(jet_pt);

    // Putting in smeared and/or scaled objects will cause that to be reflected in MET
    m_systUtil->setElectronParameters(el_smeared_pt, el_eta, el_phi,
            el_MET_wet, el_MET_wpx, el_MET_wpy,
            el_MET_statusWord);
    m_systUtil->setPhotonParameters(ph_smeared_pt, ph_eta, ph_phi,
            ph_MET_wet, ph_MET_wpx, ph_MET_wpy,
            ph_MET_statusWord);
    m_systUtil->setTauParameters(tau_pt, tau_eta, tau_phi,
            tau_MET_wet, tau_MET_wpx, tau_MET_wpy,
            tau_MET_statusWord);
    m_systUtil->setMuonParameters(mu_smeared_pt, mu_eta, mu_phi,
            mu_MET_wet, mu_MET_wpx, mu_MET_wpy,
            mu_MET_statusWord);
    // In this instance there is an overloaded version of setExtraMuonParameters that accepts smeared pTs for spectro
    m_systUtil->setExtraMuonParameters(mu_smeared_ms_pt, mu_ms_theta, mu_ms_phi);

    m_systUtil->setMETTerm(METUtil::RefTau, MET_RefTau_etx, MET_RefTau_ety, MET_RefTau_sumet);
    m_systUtil->setMETTerm(METUtil::RefMuon, MET_RefMuon_etx, MET_RefMuon_ety, MET_RefMuon_sumet);
    m_systUtil->setMETTerm(METUtil::SoftJets, MET_SoftJets_etx, MET_SoftJets_ety, MET_SoftJets_sumet);
    //   m_systUtil->setMETTerm(METUtil::CellOut, MET_CellOut_etx, MET_CellOut_ety, MET_CellOut_sumet);
    m_systUtil->setMETTerm(METUtil::CellOutEflow, MET_CellOut_Eflow_etx, MET_CellOut_Eflow_ety, MET_CellOut_Eflow_sumet);

    // These set up the systematic "SoftTerms_ptHard"
    m_systUtil->setNvtx(nvtxsoftmet);
    //  if(!isData)
    //    m_systUtil->setMETTerm(METUtil::Truth, MET_Truth_NonInt_etx, MET_Truth_NonInt_ety, MET_Truth_NonInt_sumet);

    m_systUtil->setObjectEnergyUncertainties(METUtil::Jets, jesUp, jesDown);
    m_systUtil->setObjectResolutionShift(METUtil::Jets, jerUp, jerDown);

    m_systUtil->setObjectEnergyUncertainties(METUtil::Electrons, eesUp, eesDown);
    m_systUtil->setObjectResolutionShift(METUtil::Electrons, eerUp, eerDown);

    m_systUtil->setObjectEnergyUncertainties(METUtil::Photons, pesUp, pesDown);
    m_systUtil->setObjectResolutionShift(METUtil::Photons, perUp, perDown);

    // Muons are complicated, and MET makes use of track, spectro, and combined quantites,
    // so we need all of their resolutions.
    // comboms reflects that it is the combined muon res affected by shifting ms res up and down.
    // comboid is for shifting the id res up and down
    m_systUtil->setObjectResolutionShift(METUtil::MuonsComboMS, cb_mermsUp, cb_mermsDown);
    m_systUtil->setObjectResolutionShift(METUtil::MuonsComboID, cb_meridUp, cb_meridDown);
    m_systUtil->setObjectResolutionShift(METUtil::SpectroMuons, mermsUp, mermsDown);

    // For now the mes only affects combined
    m_systUtil->setObjectEnergyUncertainties(METUtil::Muons, mesUp, mesDown);

    // Taus are just in as an example
    m_systUtil->setObjectEnergyUncertainties(METUtil::Taus, tesUp, tesDown);

    // Fill histograms

    METObject met_RefFinal = m_systUtil->getMissingET(METUtil::RefFinal);
    h_RefFinal->Fill(met_RefFinal.et()/1000);
    // Jet systematics
    METObject met_RefFinal_JESUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::JESUp);
    h_RefFinal_JESUp->Fill(met_RefFinal_JESUp.et()/1000);
    METObject met_RefFinal_JESDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::JESDown);
    h_RefFinal_JESDown->Fill(met_RefFinal_JESDown.et()/1000);
    METObject met_RefFinal_JERUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::JERUp);
    h_RefFinal_JERUp->Fill(met_RefFinal_JERUp.et()/1000);
    // Electron systematics
    METObject met_RefFinal_EESUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::EESUp);
    h_RefFinal_EESUp->Fill(met_RefFinal_EESUp.et()/1000);
    METObject met_RefFinal_EESDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::EESDown);
    h_RefFinal_EESDown->Fill(met_RefFinal_EESDown.et()/1000);
    METObject met_RefFinal_EERUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::EERUp);
    h_RefFinal_EERUp->Fill(met_RefFinal_EERUp.et()/1000);
    METObject met_RefFinal_EERDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::EERDown);
    h_RefFinal_EERDown->Fill(met_RefFinal_EERDown.et()/1000);
    // Photon systematics
    METObject met_RefFinal_PESUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::PESUp);
    h_RefFinal_PESUp->Fill(met_RefFinal_PESUp.et()/1000);
    METObject met_RefFinal_PESDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::PESDown);
    h_RefFinal_PESDown->Fill(met_RefFinal_PESDown.et()/1000);
    METObject met_RefFinal_PERUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::PERUp);
    h_RefFinal_PERUp->Fill(met_RefFinal_PERUp.et()/1000);
    METObject met_RefFinal_PERDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::PERDown);
    h_RefFinal_PERDown->Fill(met_RefFinal_PERDown.et()/1000);
    // Muon systematics
    METObject met_RefFinal_MESUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::MESUp);
    h_RefFinal_MESUp->Fill(met_RefFinal_MESUp.et()/1000);
    METObject met_RefFinal_MESDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::MESDown);
    h_RefFinal_MESDown->Fill(met_RefFinal_MESDown.et()/1000);
    METObject met_RefFinal_MERIDUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::MERIDUp);
    h_RefFinal_MERIDUp->Fill(met_RefFinal_MERIDUp.et()/1000);
    METObject met_RefFinal_MERIDDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::MERIDDown);
    h_RefFinal_MERIDDown->Fill(met_RefFinal_MERIDDown.et()/1000);
    METObject met_RefFinal_MERMSUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::MERMSUp);
    h_RefFinal_MERMSUp->Fill(met_RefFinal_MERMSUp.et()/1000);
    METObject met_RefFinal_MERMSDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::MERMSDown);
    h_RefFinal_MERMSDown->Fill(met_RefFinal_MERMSDown.et()/1000);

    // ResoSoftTerms uses gRandom for smearing. Set the seed here however you like.
    if(isData) gRandom->SetSeed(UInt_t(RunNumber * EventNumber));
    else gRandom->SetSeed(UInt_t(mc_channel_number * EventNumber));
    // Soft terms systematics
    METObject met_RefFinal_ScaleSoftTermsUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsUp);
    h_RefFinal_ScaleSoftTermsUp->Fill(met_RefFinal_ScaleSoftTermsUp.et()/1000);
    METObject met_RefFinal_ScaleSoftTermsDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsDown);
    h_RefFinal_ScaleSoftTermsDown->Fill(met_RefFinal_ScaleSoftTermsDown.et()/1000);
    METObject met_RefFinal_ResoSoftTermsUp = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUp);
    h_RefFinal_ResoSoftTermsUp->Fill(met_RefFinal_ResoSoftTermsUp.et()/1000);
    METObject met_RefFinal_ResoSoftTermsDown = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsDown);
    h_RefFinal_ResoSoftTermsDown->Fill(met_RefFinal_ResoSoftTermsDown.et()/1000);
    METObject met_RefFinal_ScaleSoftTermsUp_ptHard = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsUp_ptHard);
    // Alternate soft terms systematics
    h_RefFinal_ScaleSoftTermsUp_ptHard->Fill(met_RefFinal_ScaleSoftTermsUp_ptHard.et()/1000);
    METObject met_RefFinal_ScaleSoftTermsDown_ptHard = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ScaleSoftTermsDown_ptHard);
    h_RefFinal_ScaleSoftTermsDown_ptHard->Fill(met_RefFinal_ScaleSoftTermsDown_ptHard.et()/1000);
    METObject met_RefFinal_ResoSoftTermsUp_ptHard = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUp_ptHard);
    h_RefFinal_ResoSoftTermsUp_ptHard->Fill(met_RefFinal_ResoSoftTermsUp_ptHard.et()/1000);
    METObject met_RefFinal_ResoSoftTermsUpDown_ptHard = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsUpDown_ptHard);
    h_RefFinal_ResoSoftTermsUpDown_ptHard->Fill(met_RefFinal_ResoSoftTermsUpDown_ptHard.et()/1000);
    METObject met_RefFinal_ResoSoftTermsDownUp_ptHard = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsDownUp_ptHard);
    h_RefFinal_ResoSoftTermsDownUp_ptHard->Fill(met_RefFinal_ResoSoftTermsDownUp_ptHard.et()/1000);
    METObject met_RefFinal_ResoSoftTermsDown_ptHard = m_systUtil->getMissingET(METUtil::RefFinal,METUtil::ResoSoftTermsDown_ptHard);
    h_RefFinal_ResoSoftTermsDown_ptHard->Fill(met_RefFinal_ResoSoftTermsDown_ptHard.et()/1000);

    // Print out some information on each event
    if(verbose) {
        cout << "Demonstration of smearing and systematics" << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "All these are the scalar MET term" << endl;
        cout << "RefEle (smeared): " << m_systUtil->getMissingET(METUtil::RefEle).et() << endl;
        //     cout << "RefGamma: " << m_systUtil->getMissingET(METUtil::RefGamma).et() << endl;
        cout << "RefTau: " << m_systUtil->getMissingET(METUtil::RefTau).et() << endl;
        cout << "RefJet: " << m_systUtil->getMissingET(METUtil::RefJet).et() << endl;
        cout << "SoftJets: " << m_systUtil->getMissingET(METUtil::SoftJets).et() << endl;
        cout << "RefMuon: " << m_systUtil->getMissingET(METUtil::RefMuon).et() << endl;
        cout << "MuonBoy (smeared and scaled): " << m_systUtil->getMissingET(METUtil::MuonTotal).et() << endl;
        cout << "CellOut_Eflow: " << m_systUtil->getMissingET(METUtil::CellOutEflow).et() << endl;
        cout << "RefFinal: " << m_systUtil->getMissingET(METUtil::RefFinal).et() << endl;
        cout << "Truth: " << m_systUtil->getMissingET(METUtil::Truth).et() << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "HardTerms: " << m_systUtil->getMissingET(METUtil::HardTerms).et() << endl;
        cout << "HardTerms stands for the sum of Ref* and MuonBoy." << endl;
        cout << "SoftTerms: " << m_systUtil->getMissingET(METUtil::SoftTerms).et() << endl;
        cout << "SoftTerms stands for the sum of CellOut(_Eflow) and SoftJets." << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << endl;

        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "RefJet JESUp: " << m_systUtil->getMissingET(METUtil::RefJet,METUtil::JESUp).et()
            << ",  JESDown: " << m_systUtil->getMissingET(METUtil::RefJet,METUtil::JESDown).et() << endl;
        cout << "RefJet JES Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefJet,METUtil::JES).et() << endl;
        cout << "RefFinal JESUp: " << met_RefFinal_JESUp.et()
            << ", JESDown: " << met_RefFinal_JESDown.et() << endl;
        cout << "RefFinal JES Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::JES).et() << endl;
        cout << "RefJet JERUp: " << m_systUtil->getMissingET(METUtil::RefJet,METUtil::JERUp).et()
            << ", JERDown: " << m_systUtil->getMissingET(METUtil::RefJet,METUtil::JERDown).et() << endl;
        cout << "RefJet JER Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefJet,METUtil::JER).et() << endl;
        cout << "RefFinal JERUp: " << met_RefFinal_JERUp.et()
            << ", JERDown: " << m_systUtil->getMissingET(METUtil::RefFinal,METUtil::JERDown).et() << endl;
        cout << "RefFinal JER Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::JER).et() << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "RefEle EESUp: " << m_systUtil->getMissingET(METUtil::RefEle,METUtil::EESUp).et()
            << ",  EESDown: " << m_systUtil->getMissingET(METUtil::RefEle,METUtil::EESDown).et() << endl;
        cout << "RefFinal EESUp: " << met_RefFinal_EESUp.et()
            << ",  EESDown: " << met_RefFinal_EESDown.et() << endl;
        cout << "RefFinal EES Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::EES).et() << endl;
        cout << "RefEle EERUp: " << m_systUtil->getMissingET(METUtil::RefEle,METUtil::EERUp).et()
            << ",  EERDown: " << m_systUtil->getMissingET(METUtil::RefEle,METUtil::EERDown).et() << endl;
        cout << "RefFinal EERUp: " << met_RefFinal_EERUp.et()
            << ",  EERDown: " << met_RefFinal_EERDown.et() << endl;
        cout << "RefFinal EER Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::EER).et() << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "MuonBoy MESUp: " << m_systUtil->getMissingET(METUtil::MuonTotal,METUtil::MESUp).et()
            << ",  MESDown: " << m_systUtil->getMissingET(METUtil::MuonTotal,METUtil::MESDown).et() << endl;
        cout << "RefFinal MESUp: " << met_RefFinal_MESUp.et()
            << ",  MESDown: " << met_RefFinal_MESDown.et() << endl;
        cout << "RefFinal MES Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::MES).et() << endl;
        cout << "MuonBoy MERIDUp: " << m_systUtil->getMissingET(METUtil::MuonTotal,METUtil::MERIDUp).et()
            << ",  MERIDDown: " << m_systUtil->getMissingET(METUtil::MuonTotal,METUtil::MERIDDown).et() << endl;
        cout << "RefFinal MERIDUp: " << met_RefFinal_MERIDUp.et()
            << ",  MERIDDown: " << met_RefFinal_MERIDDown.et() << endl;
        cout << "RefFinal MERID Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::MERID).et() << endl;
        cout << "MuonBoy MERMSUp: " << m_systUtil->getMissingET(METUtil::MuonTotal,METUtil::MERMSUp).et()
            << ",  MERMSDown: " << m_systUtil->getMissingET(METUtil::MuonTotal,METUtil::MERMSDown).et() << endl;
        cout << "RefFinal MERMSUp: " << met_RefFinal_MERMSUp.et()
            << ",  MERMSDown: " << met_RefFinal_MERMSDown.et() << endl;
        cout << "RefFinal MERMS Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::MERMS).et() << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "RefTau TESUp: " << m_systUtil->getMissingET(METUtil::RefTau,METUtil::TESUp).et()
            << ",  TESDown: " << m_systUtil->getMissingET(METUtil::RefTau,METUtil::TESDown).et() << endl;
        cout << "RefFinal TESUp: " << m_systUtil->getMissingET(METUtil::RefFinal,METUtil::TESUp).et()
            << ",  TESDown: " << m_systUtil->getMissingET(METUtil::RefFinal,METUtil::TESDown).et() << endl;
        cout << "RefFinal TES Diff (up - Down)/none : " << m_systUtil->getMissingETDiff(METUtil::RefFinal,METUtil::TES).et() << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;


        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "Now follow the soft terms. Apart from AllClusters, these already include pileup contributions." << endl;
        cout << endl;
        cout << "AllClusters is the PLHC systematic on CellOut(_Eflow) and SoftJets." << endl;
        cout << "RefFinal AllClusters Up: " << m_systUtil->getMissingET(METUtil::RefFinal,METUtil::AllClustersUp).et()
            << ", RefFinal AllClusters Down: " << m_systUtil->getMissingET(METUtil::RefFinal,METUtil::AllClustersDown).et() << endl;

        cout << endl;
        cout << "These are the April 2012 systematics. For information, please see:" << endl;
        cout << "https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=161247" << endl;
        cout << endl;

        cout << "ScaleSoftTerms is the systematic on the scale CellOut(_Eflow) and SoftJets." << endl;

        // ResoSoftTerms uses gRandom for smearing. Set the seed here however you like.
        if(isData) gRandom->SetSeed(UInt_t(RunNumber * EventNumber));
        else gRandom->SetSeed(UInt_t(mc_channel_number * EventNumber));

        cout << "RefFinal ScaleSoftTerms Up: " << met_RefFinal_ScaleSoftTermsUp.et()
            << ", RefFinal ScaleSoftTerms Down: " << met_RefFinal_ScaleSoftTermsDown.et() << endl;
        cout << "ResoSoftTerms is the systematic on the scale CellOut(_Eflow) and SoftJets." << endl;
        cout << "RefFinal ResoSoftTerms Up: " << met_RefFinal_ResoSoftTermsUp.et()
            << ", RefFinal ResoSoftTerms Down: " << met_RefFinal_ResoSoftTermsDown.et() << endl;

        cout << "ScaleSoftTerms_ptHard is the systematic on the scale CellOut(_Eflow) and SoftJets." << endl;
        cout << "RefFinal ScaleSoftTerms_ptHard Up: " << met_RefFinal_ScaleSoftTermsUp_ptHard.et()
            << ", RefFinal ScaleSoftTerms_ptHard Down: " << met_RefFinal_ScaleSoftTermsDown_ptHard.et() << endl;
        cout << "ResoSoftTerms is the systematic on the scale CellOut(_Eflow) and SoftJets." << endl;
        cout << "This is parameterised in terms of longitudinal and transverse components, which can be varied coherently or anti-coherently." << endl;
        cout << "RefFinal ResoSoftTerms_ptHard Up: " << met_RefFinal_ResoSoftTermsUp_ptHard.et()
            << ", RefFinal ResoSoftTerms_ptHard Down: " << met_RefFinal_ResoSoftTermsDown_ptHard.et() << endl;
        cout << "RefFinal ResoSoftTerms_ptHard UpDown: " << met_RefFinal_ResoSoftTermsUpDown_ptHard.et()
            << ", RefFinal ResoSoftTerms_ptHard DownUp: " << met_RefFinal_ResoSoftTermsDownUp_ptHard.et() << endl;

        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << "Combined errors, giving an uncertainty on MET" << endl;
        METObject smearedMET = m_systUtil->getMissingET(METUtil::RefFinal);
        cout << "RefFinal MET = " << smearedMET.et() << " +- " << m_systUtil->absDeltaMissingET(METUtil::RefFinal).et()
            << " (" << 100*m_systUtil->deltaMissingET(METUtil::RefFinal).et() << "%)" << endl;
        cout << "+++++++++++++++++++++++++++++" << endl;
        cout << endl;
    }

}

void Example::Terminate()
{
    m_outfile->Write();
    m_outfile->Close();
}

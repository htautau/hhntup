from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
import os

from ROOT import std, TLorentzVector, TFile, Double, Long
from higgstautau.mixins import *
from externaltools import MissingETUtility
from ROOT import METUtility
from ROOT import METUtil


#################################################
# Muon Pt Smearing
#################################################

class MuonPtSmearing(EventFilter):
    """
    Smears the Muon Pt using the official ATLAS tool
    """

    def __init__(self, datatype, year, tool, **kwargs):

        super(MuonPtSmearing, self).__init__(**kwargs)
        self.datatype = datatype

        self.year = year
        self.tool = tool
        
    def passes(self, event):

        ## No smearing applied on data muons
        if self.datatype == datasets.DATA: return True

        ## Select correct smearing tool configuration
        
        for i, mu in enumerate(event.muons):

            ## Obtain parameters for correction
            charge = mu.charge
            eta    = mu.fourvect.Eta()
            pt     = mu.fourvect.Pt()
            pt_ms  = 0.0
            pt_id  = 0.0
            if mu.me_qoverp_exPV != 0:
                pt_ms  = sin(mu.me_theta_exPV)/abs(mu.me_qoverp_exPV)
            if mu.id_qoverp_exPV:
                pt_id  = sin(mu.id_theta_exPV)/abs(mu.id_qoverp_exPV)

            pt_standalone_forMET = sin(mu.ms_theta)/abs(mu.ms_qoverp)
            pt_forMET = -1.0

            
            ## Seed with event number, reproducible smear for different analyses
            self.tool.SetSeed(event.EventNumber, i)

            ## Pass parameters, get smeared Pt
            pt_smear = -1

            if mu.isCombinedMuon:

                self.tool.Event(pt_ms, pt_id, pt, eta, charge)
                pt_smear = self.tool.pTCB()
                if self.year == 2012:
                    self.tool.Event(pt_standalone_forMET, eta, 'MS', charge)
                    pt_forMET = self.tool.pTMS()
                if self.year == 2011:
                    pt_forMET = pt_smear
                
            elif mu.isSegmentTaggedMuon:
                self.tool.Event(pt_id, eta, 'ID', charge)
                pt_smear = self.tool.pTID()
                if self.year == 2012:
                    self.tool.Event(pt_standalone_forMET, eta, 'MS', charge)
                    pt_forMET = self.tool.pTMS()
                if self.year == 2011:
                    pt_forMET = pt_smear

            else:
                self.tool.Event(pt_ms, eta, 'MS', charge)
                pt_smear = self.tool.pTMS()
                if self.year == 2012:
                    pt_forMET = pt_smear
                if self.year == 2011:
                    self.tool.Event(pt_standalone_forMET, eta, 'MS', charge)
                    pt_forMET = self.tool.pTMS()

            self.tool.Event(pt_standalone_forMET, eta, 'MS', charge);
            pt_standalone_forMET  = self.tool.pTMS();
                    
            ## Adjust Pt in transient D3PD
            mu.etcone30  = pt_forMET
            mu.ptcone30 = pt_standalone_forMET
            mu.pt = pt_smear

            ## Adjust Pt of the muon 4-vector
            newMufourvect = TLorentzVector()
            newMufourvect.SetPtEtaPhiM(pt_smear, eta, mu.fourvect.Phi(), mu.fourvect.M())
            setattr(mu, 'fourvect', newMufourvect)

            # ## Adjust MET accordingly
            # px = pt*cos(mu.phi)
            # py = pt*sin(mu.phi)
            
            # if self.year == 2011:
            #     px_smear = pt_standalone_forMET*cos(mu.phi)
            #     py_smear = pt_standalone_forMET*sin(mu.phi)

            #     event.MET_RefFinal_BDTMedium_etx += (px-px_smear)*mu.MET_BDTMedium_wpx[0]
            #     event.MET_RefFinal_BDTMedium_ety += (py-py_smear)*mu.MET_BDTMedium_wpy[0]
            #     event.MET_RefFinal_BDTMedium_sumet -= (pt-pt_smear)*mu.MET_BDTMedium_wet[0]

            # if self.year == 2012:
            #     px_smear = pt_standalone_forMET*cos(mu.phi)
            #     py_smear = pt_standalone_forMET*sin(mu.phi)

            #     event.MET_RefFinal_STVF_etx += (px-px_smear)*mu.MET_STVF_wpx[0]
            #     event.MET_RefFinal_STVF_ety += (py-py_smear)*mu.MET_STVF_wpy[0]
            #     event.MET_RefFinal_STVF_sumet -= (pt-pt_smear)*mu.MET_STVF_wet[0]

        return True


#################################################
# Muon isolation correction
#################################################
class MuonIsoCorrection(EventFilter):
    """
    Corrects the muon isolation in 2012
    """

    def __init__(self, datatype, year, tool, **kwargs):

        super(MuonIsoCorrection, self).__init__(**kwargs)
        self.datatype = datatype
        self.year = year
        self.tool = tool

    def passes(self, event):

        ## Apply only in 2012
        if self.year < 2012: return True

        ##Calculate HSG3 nvtx
        hsg3_nvtx = 0
        
        for vx in event.vertices:
            if vx.nTracks >= 3:
                hsg3_nvtx += 1
            
        ## Apply the correction
        for mu in event.muons:

            newEtcone20 = self.tool.CorrectEtCone(mu.etcone20, float(hsg3_nvtx), mu.fourvect.Eta(), 'cone20Comb2012')
            mu.etcone20 = newEtcone20

        return True



#################################################
# Correct MET
#################################################
class METCorrection(EventFilter):
    """
    Adjusts the MET according to new adjusted Pt of objects
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Missing_ET_Correction_due_to_Pt
    """

    def __init__(self, datatype, year, **kwargs):

        super(METCorrection, self).__init__(**kwargs)
        self.tool = METUtility()
        self.datatype = datatype
        self.year = year

    def passes(self, event):

        if self.datatype == datasets.DATA: return True

        ## Collect jet Info
        jets_calib_pt=event.jet_pt
        jets_calib_eta=event.jet_eta
        jets_calib_phi=event.jet_phi
        jets_old_E=event.jet_PhiOriginEM
        jets_old_pt=event.jet_EtaOriginEM
            
        if self.year == 2011:
            jets_MET_wpx=event.jet_AntiKt4LCTopo_MET_BDTMedium_wpx
            jets_MET_wpy=event.jet_AntiKt4LCTopo_MET_BDTMedium_wpy
            jets_MET_wet=event.jet_AntiKt4LCTopo_MET_BDTMedium_wet
            jets_MET_status=event.jet_AntiKt4LCTopo_MET_BDTMedium_statusWord
        if self.year == 2012:
            jets_MET_wpx=event.jet_AntiKt4LCTopo_MET_STVF_wpx
            jets_MET_wpy=event.jet_AntiKt4LCTopo_MET_STVF_wpy
            jets_MET_wet=event.jet_AntiKt4LCTopo_MET_STVF_wet
            jets_MET_status=event.jet_AntiKt4LCTopo_MET_STVF_statusWord
                
        self.tool.setObjects(METUtil.Jets,
                             jets_calib_pt,
                             jets_calib_eta,
                             jets_calib_phi,
                             jets_old_E,
                             jets_MET_wet,
                             jets_MET_wpx, 
                             jets_MET_wpy,
                             jets_MET_status)

        self.tool.setOriJetParameters(jets_old_pt)


        
        ## Collect muon Info
        muons_calib_pt=event.mu_staco_etcone30
        muons_calib_eta=event.mu_staco_eta
        muons_calib_phi=event.mu_staco_phi
        muons_staco_charge=event.mu_staco_charge

        muons_ms_pt=event.mu_staco_ptcone30
        muons_ms_theta=event.mu_staco_ms_theta
        muons_ms_phi=event.mu_staco_ms_phi

        if self.year == 2011:
            muons_MET_wet=event.mu_staco_MET_BDTMedium_wet
            muons_MET_wpx=event.mu_staco_MET_BDTMedium_wpx
            muons_MET_wpy=event.mu_staco_MET_BDTMedium_wpy
            muons_MET_status=event.mu_staco_MET_BDTMedium_statusWord

        if self.year == 2012:
            muons_MET_wet=event.mu_staco_MET_STVF_wet
            muons_MET_wpx=event.mu_staco_MET_STVF_wpx
            muons_MET_wpy=event.mu_staco_MET_STVF_wpy
            muons_MET_status=event.mu_staco_MET_STVF_statusWord
        
        self.tool.setMuonParameters(muons_calib_pt,
                                    muons_calib_eta,
                                    muons_calib_phi,
                                    muons_MET_wet,
                                    muons_MET_wpx,
                                    muons_MET_wpy,
                                    muons_MET_status)

        self.tool.setExtraMuonParameters(muons_ms_pt,
                                         muons_ms_theta,
                                         muons_ms_phi)

        ## Collect electron Info
        electrons_calib_pt=event.el_cl_pt
        electrons_calib_eta=event.el_eta
        electrons_calib_phi=event.el_phi

        if self.year == 2011:
            electrons_MET_wet=event.el_MET_BDTMedium_wet
            electrons_MET_wpx=event.el_MET_BDTMedium_wpx
            electrons_MET_wpy=event.el_MET_BDTMedium_wpy
            electrons_MET_status=event.el_MET_BDTMedium_statusWord

        if self.year == 2012:
            electrons_MET_wet=event.el_MET_STVF_wet
            electrons_MET_wpx=event.el_MET_STVF_wpx
            electrons_MET_wpy=event.el_MET_STVF_wpy
            electrons_MET_status=event.el_MET_STVF_statusWord

        self.tool.setElectronParameters(electrons_calib_pt,
                                        electrons_calib_eta,
                                        electrons_calib_phi,
                                        electrons_MET_wet,
                                        electrons_MET_wpx,
                                        electrons_MET_wpy,
                                        electrons_MET_status)

        ## Collect electron Info
        taus_calib_pt=event.tau_pt
        taus_calib_eta=event.tau_eta
        taus_calib_phi=event.tau_phi

        if self.year == 2011:
            taus_MET_wet=event.tau_MET_BDTMedium_wet
            taus_MET_wpx=event.tau_MET_BDTMedium_wpx
            taus_MET_wpy=event.tau_MET_BDTMedium_wpy
            taus_MET_status=event.tau_MET_BDTMedium_statusWord

        if self.year == 2012:
            taus_MET_wet=event.tau_MET_STVF_wet
            taus_MET_wpx=event.tau_MET_STVF_wpx
            taus_MET_wpy=event.tau_MET_STVF_wpy
            taus_MET_status=event.tau_MET_STVF_statusWord

        self.tool.setTauParameters(taus_calib_pt,
                                   taus_calib_eta,
                                   taus_calib_phi,
                                   taus_MET_wet,
                                   taus_MET_wpx,
                                   taus_MET_wpy,
                                   taus_MET_status)

        ## Collect MET Terms
        if self.year == 2011:
            self.tool.setMETTerm(METUtil.RefGamma,
                                 event.MET_RefGamma_BDTMedium_etx,
                                 event.MET_RefGamma_BDTMedium_ety,
                                 event.MET_RefGamma_BDTMedium_sumet)
            
            self.tool.setMETTerm(METUtil.CellOutEflow,
                                 event.MET_CellOut_BDTMedium_etx,
                                 event.MET_CellOut_BDTMedium_ety,
                                 event.MET_CellOut_BDTMedium_sumet)
            
            self.tool.setMETTerm(METUtil.SoftJets,
                                 event.MET_SoftJets_BDTMedium_etx,
                                 event.MET_SoftJets_BDTMedium_ety,
                                 event.MET_SoftJets_BDTMedium_sumet)
            
            self.tool.setMETTerm(METUtil.RefMuon,
                                 event.MET_RefMuon_Staco_BDTMedium_etx,
                                 event.MET_RefMuon_Staco_BDTMedium_ety,
                                 event.MET_RefMuon_Staco_BDTMedium_sumet)

        if self.year == 2012:
            self.tool.setMETTerm(METUtil.RefGamma,
                                 event.MET_RefGamma_STVF_etx,
                                 event.MET_RefGamma_STVF_ety,
                                 event.MET_RefGamma_STVF_sumet)
            
            self.tool.setMETTerm(METUtil.CellOutEflow,
                                 event.MET_CellOut_STVF_etx,
                                 event.MET_CellOut_STVF_ety,
                                 event.MET_CellOut_STVF_sumet)
            
            self.tool.setMETTerm(METUtil.RefMuon,
                                 event.MET_RefMuon_Staco_STVF_etx,
                                 event.MET_RefMuon_Staco_STVF_ety,
                                 event.MET_RefMUon_Staco_STVF_sumet)

        newMET = self.tool.getMissingET(METUtil.RefFinal)
            
        if self.year == 2011:
            event.MET_RefFinal_BDTMedium_etx = newMET.etx()
            event.MET_RefFinal_BDTMedium_ety = newMET.ety()
            event.MET_RefFinal_BDTMedium_sumet = newMET.sumet()

        if self.year == 2012:
            event.MET_RefFinal_STVF_etx = newMET.etx()
            event.MET_RefFinal_STVF_ety = newMET.ety()
            event.MET_RefFinal_STVF_sumet = newMET.sumet()

        return True
                             
        
    
        


#################################################
# Electron Energy Rescaling
#################################################

class EgammaERescaling(EventFilter):
    """
    Rescales/smeares the electron energy using the official ATLAS tool.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyScaleResolutionRecommendations
    """

    def __init__(self, datatype, year, tool, **kwargs):

        super(EgammaERescaling, self).__init__(**kwargs)
        self.datatype = datatype
        self.year = year
        self.tool = tool
        if self.year == 2011:
            self.tool.useDefaultCalibConstants('2011')
        if self.year == 2012:
            self.tool.useDefaultCalibConstants('2012')

    def passes(self, event):

        for el in event.electrons:

            raw_e  = el.cl_E
            raw_et = getattr(el,'fourvect').Pt()
            corrected_e  = raw_e
            corrected_et = raw_et
            cl_eta = el.cl_eta
            trk_eta = el.tracketa
            cl_phi = el.cl_phi
            
            ## Seed with event number, reproducible smear for different analyses
            self.tool.SetRandomSeed(abs(int(5*event.RunNumber+event.EventNumber+(cl_phi+pi)*1000000.)))

            if self.year == 2011:
                ## Calibration for electrons in the transition region in Data and MC
                scaleFactorForTransitionRegion = self.tool.applyMCCalibrationMeV(cl_eta, raw_et, 'ELECTRON')
                corrected_e  *= scaleFactorForTransitionRegion
                corrected_et = corrected_e/cosh(trk_eta)

                if self.datatype == datasets.MC or self.datatype == datasets.EMBED:
                    ## Smearing correction in MC
                    sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                    smearFactor = self.tool.getSmearingCorrectionMeV(cl_eta, corrected_e, sys, False, '2011')
                    corrected_e  *= smearFactor
                    corrected_et = corrected_e/cosh(trk_eta)
                else:
                    ## Calibration correction in Data
                    sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                    scaleFactor = 1
                    if corrected_e > 0:
                        scaleFactor = self.tool.applyEnergyCorrectionMeV(cl_eta, cl_phi, corrected_e, corrected_et, sys, 'ELECTRON') / corrected_e
                    corrected_e  *= scaleFactor
                    corrected_et = corrected_e/cosh(trk_eta)

            if self.year == 2012:
                if self.datatype == datasets.DATA:
                    corrected_e = self.tool.applyEnergyCorrectionMeV(cl_eta, cl_phi, raw_e, raw_et, 0, 'ELECTRON')
                    corrected_et = corrected_e/cosh(trk_eta)
                else:
                    corrected_e = self.tool.applyEnergyCorrectionMeV(cl_eta, cl_phi, raw_e, raw_et, 0, 'ELECTRON')
                    corrected_e *= self.tool.getSmearingCorrectionMeV(cl_eta, corrected_e, 0, False)
                    corrected_et = corrected_e/cosh(trk_eta)


            ## Modify E and Et in transient D3PD
            el.cl_E  = corrected_e
            el.cl_et = corrected_et
            el.cl_pt = corrected_et

            ## Modify the fourvector which is used for the rest of the tree filling
            el_eta = el.fourvect.Eta()
            el_phi = el.fourvect.Phi()
            newElfourvect = TLorentzVector()
            newElfourvect.SetPtEtaPhiE(corrected_et, el_eta, el_phi, corrected_e)
            setattr(el, 'fourvect', newElfourvect)
            

            # ## Adjust MET accordingly
            # px = raw_et * cos(el.phi)
            # py = raw_et * sin(el.phi)
            # corrected_px = corrected_et * cos(el.phi)
            # corrected_py = corrected_et * sin(el.phi)

            # if self.year == 2011:
            #     event.MET_RefFinal_BDTMedium_etx += ( px - corrected_px ) * el.MET_BDTMedium_wpx[0]
            #     event.MET_RefFinal_BDTMedium_ety += ( py - corrected_py ) * el.MET_BDTMedium_wpy[0]
            #     event.MET_RefFinal_BDTMedium_sumet -= ( raw_et - corrected_et ) * el.MET_BDTMedium_wet[0]
            # if self.year == 2012:
            #     event.MET_RefFinal_STVF_etx += ( px - corrected_px ) * el.MET_STVF_wpx[0]
            #     event.MET_RefFinal_STVF_ety += ( py - corrected_py ) * el.MET_STVF_wpy[0]
            #     event.MET_RefFinal_STVF_sumet -= ( raw_et - corrected_et ) * el.MET_STVF_wet[0]
                

        return True


#################################################
# Electron Isolation correction
#################################################
    
class ElectronIsoCorrection(EventFilter):
    """
    Correction for electron calorimeter isolation variables
    """
    
    def __init__(self, datatype, year, tool, **kwargs):

        self.datatype = datatype
        self.year = year
        self.tool = tool  
        super(ElectronIsoCorrection, self).__init__(**kwargs)

    def passes(self, event):

        nPV = 0
        for vxp in event.vertices:
            if vxp.nTracks >= 2: nPV += 1
        
        for el in event.electrons:
            E = el.fourvect.E()
            EtaS2 = el.etas2
            EtaP  = el.etap
            cl_eta = el.cl_eta
            EtCone20 = el.Etcone20
            newEtCone20 = EtCone20
            
            if self.year == 2011:
                if self.datatype == datasets.MC:
                    newElectronEtCone20 = self.tool.GetPtNPVCorrectedIsolation(nPV,
                                                                               E,
                                                                               EtaS2,
                                                                               EtaP,
                                                                               cl_eta,
                                                                               20,
                                                                               True,
                                                                               EtCone20)
                else:
                    newElectronEtCone20 = self.tool.GetPtNPVCorrectedIsolation(nPV,
                                                                               E,
                                                                               EtaS2,
                                                                               EtaP,
                                                                               cl_eta,
                                                                               20,
                                                                               False,
                                                                               EtCone20)
                el.Etcone20 = newEtCone20

            if self.year == 2012:
                newTopoEtCone20 = self.tool.GetPtEDCorrectedTopoIsolation(el.ED_median,
                                                                          E,
                                                                          EtaS2,
                                                                          EtaP,
                                                                          cl_eta,
                                                                          20,
                                                                          self.datatype == datasets.MC,
                                                                          el.topoEtcone20,
                                                                          False,
                                                                          CaloIsoCorrection.ELECTRON,
                                                                          CaloIsoCorrection.REL17)
                
                el.topoEtcone20 = newTopoEtCone20
                
        return True


#################################################
# Muon isolation efficiency
#################################################
def MuonIsoEffCorrection(tool, event, datatype, year):
    """
    Apply a correction on the muon isolation, returns nominal, -1sigma, +1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    if datatype == datasets.DATA or year < 2012: return 1.0, 1.0, 1.0

    mu = event.muons[0]

    ##Calculate HSG3 nvtx
    hsg3_nvtx = 0
        
    for vx in event.vertices:
        if vx.nTracks >= 3:
            hsg3_nvtx += 1

    b = tool.FindBin(mu.fourvect.Pt()/1000, hsg3_nvtx)
    weight = tool.GetBinContent(b)

    return weight, 0.98, 1.02
    

    


#################################################
# Tau/Electron Misidentification correction
#################################################
def TauEfficiencySF(event, datatype, year):
    """
    Apply Tau Efficiency Scale Factor correction
    Returns, nominal, -1sigma, +1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    sf, errup, errdown = 1.0, 1.0, 1.0
        
    # Apply only on MC
    if datatype != datasets.MC: return sf, errup, errdown

    errup, errdown = 0.0, 0.0

    #Correct 1p taus/ electron fake rate
    #https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012#Systematics_for_Electron_Misiden

    tau = event.taus[0]
    
    eta = abs(tau.eta)
    
    if tau.numTrack == 1 and tau.fourvect.Pt() > 20*GeV:
        nMC = event.mc_n
        for i in range(0, nMC):
            if abs(event.mc_pdgId[i]) == 11 and event.mc_pt[i] > 8*GeV and event.mc_status[i] == 1:
                if utils.dR(event.mc_eta[i], event.mc_phi[i], tau.eta, tau.phi) > 0.2: continue
                
                if year == 2011:
                    if eta < 1.37:
                        sf, errup, errdown = 1.64, 0.81, 0.81 
                    elif eta < 1.52:
                        sf, errup, errdown = 1.00, 1.00, 1.00
                    elif eta < 2.0:
                        sf, errup, errdown = 0.71, 0.63, 0.63
                    else:
                        sf, errup, errdown = 2.90, 1.42, 1.42
                        
                if year == 2012:
                    if eta < 0.05:
                        sf, errup, errdown = 0.856, 0.154, 0.154
                    elif eta < 1.37:
                        sf, errup, errdown = 1.119, 0.200, 0.200
                    elif eta < 2.00:
                        sf, errup, errdown = 1.402, 0.307, 0.307
                    elif eta < 2.30:
                        sf, errup, errdown = 1.579, 1.176, 1.176
                    elif eta < 2.47:
                        sf, errup, errdown = 2.669, 2.669, 2.936
                    else:
                        sf, errup, errdown = 21.688, 21.688, 22.277

    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown



#################################################
# Tau ID correction
#################################################
def TauIDSF(event, datatype, year):
    """
    Applies the tau ID scale factor
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    ## Apply only on MC and Embedding, initialze numbers to return
    sf, errup, errdown = 1.0, 1.0, 1.0
    if datatype == datasets.DATA: return sf, errup, errdown
    errup, errdown = 0.0, 0.0

    for tau in event.taus:
        ntracks = tau.numTrack
        eta     = abs(tau.eta)

        if year == 2011:
            ## https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012
            ## https://twiki.cern.ch/twiki/pub/AtlasProtected/TauSystematicsWinterConf2012/tauWG_wtaunu2012winter_01.27.2012_soshi.pdf
            if ntracks == 1:
                sf, errup, errdown = 0.9210, 0.0421, 0.0421
            if ntracks == 3:
                sf, errup, errdown = 0.9849, 0.08, 0.08

        if year == 2012:
            ## https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsSummerConf2012
            if ntracks == 1:
                if eta < 1.5:
                    sf, errup, errdown = 0.992, 0.03, 0.03
                if eta >= 1.5:
                    sf, errup, errdown = 0.952, 0.04, 0.04
            if ntracks == 3:
                if eta < 1.5:
                    sf, errup, errdown = 1.073, 0.071, 0.071
                if eta >= 1.5:
                    sf, errup, errdown = 0.99, 0.10, 0.10

    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown


#################################################
# Muon Efficiency corrections
#################################################

## Scale factors for single muon triggers
def MuonSF(tool, tool2, event, datatype, pileup_tool, year, runNumber, isLTT, TrigMuonEff):
    """
    Apply Muon Efficiency correction
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    sf, errup, errdown = 1.0, 1.0, 1.0
    
    # Apply only on MC and Embedding
    if datatype == datasets.DATA: sf, errup, errdown

    errup, errdown = 0.0, 0.0

    muon = event.muons[0].fourvect
    muon_charge = int(event.muons[0].charge)
    muonSF = tool

    if year == 2011:
        sf = muonSF.scaleFactor(muon)
        sf_err_stat = muonSF.scaleFactorUncertainty(muon)
        sf_err_syst = muonSF.scaleFactorSystematicUncertainty(muon)
        errup = sqrt(sf_err_stat**2 + sf_err_syst**2)
        errdown = errup

        if datatype == datasets.EMBED and not isLTT:
            period = None
            if runNumber <  178044: period = TrigMuonEff.perUnDefined
            elif runNumber <= 186493: period = TrigMuonEff.per2011B_I
            elif runNumber <= 189090 or (runNumber >= 189639 and runNumber <= 191933): period = TrigMuonEff.per2011J_MwoL3_L4
            elif runNumber >= 189184 and runNumber <= 189610: period = TrigMuonEff.per2011L3_L4

            trigsf = tool2.MuEff( period, False, muon, 1)
            sf *= trigsf.first
            errup = sqrt(errup**2 + trigsf.second**2)
            errdown = errup
    
    if year == 2012:
        sf = muonSF.scaleFactor(muon_charge, muon)
        sf_err_stat = muonSF.scaleFactorUncertainty(muon_charge, muon)
        sf_err_syst = muonSF.scaleFactorSystematicUncertainty(muon_charge, muon)
        errup = sqrt(sf_err_stat**2 + sf_err_syst**2)
        errdown = errup
        
    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown

## Scale factor for lephad triggers

def MuonLTTSF(tool, event, datatype, year, pileup_tool):
    """
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    sf, errup, errdown = 1.0, 0.0, 0.0

    mu = event.muons[0].fourvect
    runNumber = event.RunNumber

    if year == 2011:
        if datatype == datasets.MC:
            sf      = tool.getSFMuon(mu, runNumber, 0)
            errup   = tool.getSFMuon(mu, runNumber, 1)
            errdown = tool.getSFMuon(mu, runNumber, -1)
        else:
            sf      = tool.getDataEffMuon(mu, runNumber, 0)
            errup   = tool.getDataEffMuon(mu, runNumber, 1)
            errdown = tool.getDataEffMuon(mu, runNumber, -1)


    elif year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(runNumber)
        
        if datatype == datasets.MC:
            sf      = tool.getSFMuon(mu, randomRunNumber, 0)
            errup   = tool.getSFMuon(mu, randomRunNumber, 1)
            errdown = tool.getSFMuon(mu, randomRunNumber, -1)
        else:
            sf      = tool.getDataEffMuon(mu, randomRunNumber, 0)
            errup   = tool.getDataEffMuon(mu, randomRunNumber, 1)
            errdown = tool.getDataEffMuon(mu, randomRunNumber, -1)

    if sf > 0.0:
        errup   = errup / sf
        errdown = errdown / sf

    return sf, errup, errdown
        
        


#################################################
# Electron Efficiency corrections
#################################################

## Scale factors for single electron triggers
def ElectronSF(tool, event, datatype, pileup_tool, year, runNumber, isLTT):

    """
    Electron efficiency correction.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EfficiencyMeasurements
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    sf, errup, errdown = 1.0, 1.0, 1.0
    
    # Apply only on MC and Embedding
    if datatype == datasets.DATA: sf, errup, errdown

    errup, errdown = 0.0, 0.0

    el = event.electrons[0]
    
    pt = el.fourvect.Pt()
    eta = el.cl_eta
    if year == 2011:
        reconstructionEffCorr = tool.scaleFactor(eta, pt, 4, 0, 6, True) # track + reco eff SF
        identificationEffCorr = tool.scaleFactor(eta, pt, 7, 0, 6, True) # ID eff SF

        if not isLTT:
            Set = None
            if datatype == datasets.EMBED:
                if runNumber < 186755: Set = 11
                elif runNumber < 187825: Set = 15
                elif runNumber < 191933: Set = 19
                    
            if datatype == datasets.MC:
                if runNumber < 186755: Set=10
                elif runNumber < 191933: Set = 18
                            
            trigEffCorr = tool.scaleFactor(eta, pt, Set, 0, 6, True)
        
    if year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(event.RunNumber)
        reconstructionEffCorr = tool.scaleFactor(eta, pt, 4, 0, 8, 1, randomRunNumber) # track + reco eff SF
        identificationEffCorr = tool.scaleFactor(eta, pt, 7, 0, 8, 1, randomRunNumber) # ID eff SF
    ## We don't use electrons in the forward region, so we don't need to apply SF for them.
    
    sf = reconstructionEffCorr.first * identificationEffCorr.first 
    if not isLTT:
        sf *= trigEffCorr.first
    reconstructionError = reconstructionEffCorr.second
    identificationError = identificationEffCorr.second
    if not isLTT:
        trigError = trigEffCorr.second
    errup = (reconstructionError**2 + identificationError**2)
    if not isLTT:
        errup += trigError**2
    errup = sqrt(errup)
    errdown = errup

    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown




#################################################
# Lepton trigger scale factors
#################################################

## Scale factor for single lepton triggers
def LeptonSLTSF(tool, event, datatype, pileup_tool, year):

    """
    Trigger efficiency correction.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Electrons
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/TrigMuonEfficiency#LeptonTriggerSF_Tool
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """
    ## From the
    ## https://twiki.cern.ch/twiki/bin/viewauth/Atlas/TrigMuonEfficiency#LeptonTriggerSF_Tool
    ## Twiki:
    ## * The RunNumber should be ideally randomly generated using the recommended pileup reweighting tool
    ##   to assure close correspondance to actual data run numbers.
    ## * The vector of muons and electrons should be filled for muons and electrons separately
    ##   using ALL good leptons that passed the object selection, not just the trigger matched ones.
    ##   Independently of this, the analysis should require that at least one lepton be trigger matched
    ##   or else a bias can be introduced (especially for single lepton analyses; the bias decreases
    ##   with the number of leptons).
    ## From the
    ## https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Electrons
    ## Twiki:
    ## * The trigger matching is not applied for this analysis because of missing information in D3PDs.
    ## From the
    ## https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Electrons
    ## Twiki:
    ## * The trigger matching is not applied for this analysis because of missing information in D3PDs,
    ##   but this is a small effect.

    sf, errup, errdown = 1.0, 1.0, 1.0
    
    # Apply only on MC and Embedding
    if datatype == datasets.DATA: return sf, errup, errdown

    errup, errdown = 0.0, 0.0
    
    v_muTLVs = std.vector(TLorentzVector)()
    v_elTLVs = std.vector(TLorentzVector)()

    for imu, mu in enumerate(event.muons):
        if imu > 0:
            print 'WARNING in LeptonSF(): using more than 1 muon in efficiency correction'
        muTLV = TLorentzVector()
        muTLV.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, mu.m)
        v_muTLVs.push_back(muTLV)

    for iel, el in enumerate(event.electrons):
        if iel > 0:
            print 'WARNING in LeptonSF(): using more than 1 electron in efficiency correction'
        elTLV = TLorentzVector()
        pt  = el.fourvect.Pt()
        eta = el.cl_eta
        phi = el.cl_phi
        E   = el.fourvect.E()
        elTLV.SetPtEtaPhiE(pt,eta,phi,E)
        v_elTLVs.push_back(elTLV)

    if year == 2011:
        randomRunNumber = event.RunNumber
    if year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(event.RunNumber)

    leptonTriggerSF = tool

    trigSF = leptonTriggerSF.GetTriggerSF(randomRunNumber, False, v_muTLVs, 1, v_elTLVs, 2, 0)
    sf = trigSF.first
    errup = trigSF.second
    errdown = errup

    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown



## Scale factor for lephad triggers
def ElectronLTTSF(tool, event, datatype, year, pileup_tool, runNumber):
    """
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    sf, errup, errdown = 1.0, 1.0, 1.0
    
    # Apply only on MC and Embedding
    if datatype == datasets.DATA: sf, errup, errdown

    errup, errdown = 0.0, 0.0

    e = event.electrons[0].fourvect
    runNumber = event.RunNumber

    if year == 2011:
        if datatype == datasets.MC:
            sf      = tool.getSFElec(e, runNumber, 0)
            errup   = tool.getSFElec(e, runNumber, 1)
            errdown = tool.getSFElec(e, runNumber, -1)
        else:
            sf      = tool.getDataEffElec(e, runNumber, 0)
            errup   = tool.getDataEffElec(e, runNumber, 1)
            errdown = tool.getDataEffElec(e, runNumber, -1)

    elif year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(runNumber)
        
        if datatype == datasets.MC:
            sf      = tool.getSFElec(e, randomRunNumber, 0)
            errup   = tool.getSFElec(e, randomRunNumber, 1)
            errdown = tool.getSFElec(e, randomRunNumber, -1)
        else:
            sf      = tool.getDataEffElec(e, randomRunNumber, 0)
            errup   = tool.getDataEffElec(e, randomRunNumber, 1)
            errdown = tool.getDataEffElec(e, randomRunNumber, -1)

    if sf > 0.0:
        errup   = errup / sf
        errdown = errdown / sf

    return sf, errup, errdown


    

#################################################
# Special LTT corrections for tau trigger
#################################################

def TauTriggerSF(tool1, tool2, tool3, Tau, nvtx, runNumber, year, lep, pileup_tool, datatype):
    
    sf, errup, errdown = 1.0, 0.0, 0.0
    tauPt  = Tau.fourvect.Pt()
    tauEta = Tau.fourvect.Eta()

    if year == 2011:
        if datatype == datasets.EMBED:
            if runNumber > 186755:
                sf    =  tool1.get3DMCEff(tauPt, tauEta , nvtx, 0)
                errup    += tool1.get3DMCEff(tauPt, tauEta , nvtx, 1)**2
                errdown  += tool1.get3DMCEff(tauPt, tauEta , nvtx, -1)**2
                sf   *= tool1.getSF(tauPt, 0)
                errup    += tool1.getSF(tauPt, 1)**2
                errdown  += tool1.getSF(tauPt, -1)**2
                errup = sqrt(errup)
                errdown = sqrt(errdown)
            else:
                sf  = tool3.getDataEff(tauPt, 0)
                errup   = tool3.getDataEff(tauPt, 1)
                errdown = tool3.getDataEff(tauPt, -1)
                
        if datatype == datasets.MC:
            if runNumber > 186755:
                sf   = tool2.getSF(tauPt, 0)
                errup    = tool2.getSF(tauPt, 1)
                errdown  = tool2.getSF(tauPt, -1)
            else:
                sf  = tool3.getSF(tauPt, 0)
                errup   = tool3.getSF(tauPt, 1)
                errdown = tool3.getSF(tauPt, -1)

    if year == 2012:
        
        randomRunNumber = pileup_tool.GetRandomRunNumber(runNumber)
        period = ''

        if 200804 < randomRunNumber <= 201556: period = 'periodA'
        if 202660 < randomRunNumber <= 209025: period = 'periodBD'
        if 209073 < randomRunNumber <= 210308: period = 'periodE'

        nprongs = ''
            
        if Tau.numTrack == 1:
            nprong = '1p'
        if Tau.numTrack == 3:
            nprong = '3p'
        
        if lep == 'mu':
            sf  = tool1.getSF(tauPt, tauEta, 0, period, nprong, 'BDTm', 'EVm')
            errup   = tool1.getSF(tauPt, tauEta, 1, period, nprong, 'BDTm', 'EVm')
            errdown = tool1.getSF(tauPt, tauEta, -1, period, nprong, 'BDTm', 'EVm')

        if lep == 'e':
            sf  = tool2.getSF(tauPt, tauEta, 0, period, nprong, 'BDTm', 'EVm')
            errup   = tool2.getSF(tauPt, tauEta, 1, period, nprong, 'BDTm', 'EVm')
            errdown = tool2.getSF(tauPt, tauEta, -1, period, nprong, 'BDTm', 'EVm')

    if sf > 0.0:
        errup   = errup / sf
        errdown = errdown / sf

    return sf, errup, errdown
            
        

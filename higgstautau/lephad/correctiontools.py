from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
import os

from externaltools import MuonMomentumCorrections
from externaltools import MuonEfficiencyCorrections
from externaltools import TrigMuonEfficiency

from externaltools import egammaAnalysisUtils

from externaltools import HSG4TriggerSF as HSG4
from externaltools import TauTriggerCorrections as TTC

from ROOT import std, TLorentzVector, TFile
from higgstautau.mixins import *


#################################################
# Muon Pt Smearing
#################################################

from ROOT import MuonSmear

muonSmear = MuonSmear.SmearingClass(
    "Data11","staco","pT","Rel17",
    MuonMomentumCorrections.RESOURCE_PATH)

muonSmear.UseScale(1)
muonSmear.UseImprovedCombine()
muonSmear.RestrictCurvatureCorrections(2.5)
muonSmear.FillScales("KC")

class MuonPtSmearing(EventFilter):
    """
    Smears the Muon Pt using the official ATLAS tool
    """

    def __init__(self, datatype, **kwargs):

        super(MuonPtSmearing, self).__init__(**kwargs)
        self.datatype = datatype

    def passes(self, event):

        if self.datatype == datasets.DATA: return True

        for mu in event.muons:

            ## Obtain parameters for correction
            charge = mu.charge
            eta    = mu.fourvect.Eta()
            pt     = mu.fourvect.Pt()
            pt_ms  = sin(mu.ms_theta)/abs(mu.ms_qoverp)
            pt_id  = sin(mu.id_theta)/abs(mu.id_qoverp)

            ## Seed with event number, reproducible smear for different analyses
            muonSmear.SetSeed(event.EventNumber)

            ## Pass parameters, get smeared Pt
            muonSmear.Event(pt_ms, pt_id, pt, eta, charge)
            pt_smear = -1

            if mu.isCombinedMuon:
                pt_smear = muonSmear.pTCB()
            else:
                pt_smear = muonSmear.pTID()

            ## Adjust Pt in transient D3PD
            mu.pt = pt_smear

            ## Adjust Pt of the muon 4-vector
            newMufourvect = TLorentzVector()
            newMufourvect.SetPtEtaPhiM(pt_smear, eta, mu.fourvect.Phi(), mu.fourvect.M())
            setattr(mu, 'fourvect', newMufourvect)

            ## Adjust MET accordingly
            px = pt*cos(mu.phi)
            py = pt*sin(mu.phi)

            px_smear = pt_smear*cos(mu.phi)
            py_smear = pt_smear*sin(mu.phi)

            event.MET_RefFinal_BDTMedium_etx += (px-px_smear)*mu.MET_BDTMedium_wpx[0]
            event.MET_RefFinal_BDTMedium_ety += (py-py_smear)*mu.MET_BDTMedium_wpy[0]
            event.MET_RefFinal_BDTMedium_sumet -= (pt-pt_smear)*mu.MET_BDTMedium_wet[0]

        return True


#################################################
# Electron Energy Rescaling
#################################################

from ROOT import eg2011

egammaER = eg2011.EnergyRescaler()
egammaER.useDefaultCalibConstants("2011")

class EgammaERescaling(EventFilter):
    """
    Rescales/smeares the electron energy using the official ATLAS tool.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyScaleResolutionRecommendations
    """

    def __init__(self, datatype, **kwargs):

        super(EgammaERescaling, self).__init__(**kwargs)
        self.datatype = datatype

    def passes(self, event):

        for el in event.electrons:

            ## Seed with event number, reproducible smear for different analyses
            egammaER.SetRandomSeed(int(5*event.RunNumber+event.EventNumber+(getattr(el,'fourvect').Phi()+pi)*1000000.))

            raw_e  = el.cl_E
            raw_et = getattr(el,'fourvect').Pt()
            corrected_e  = raw_e
            corrected_et = raw_et
            cl_eta = el.cl_eta

            ## Calibration for electrons in the transition region in Data and MC
            scaleFactorForTransitionRegion = egammaER.applyMCCalibrationMeV(cl_eta, raw_et, 'ELECTRON')
            corrected_e  *= scaleFactorForTransitionRegion
            corrected_et *= scaleFactorForTransitionRegion

            if self.datatype == datasets.MC or self.datatype == datasets.EMBED:
                ## Smearing correction in MC
                sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                smearFactor = egammaER.getSmearingCorrectionMeV(cl_eta, corrected_e, sys, False, '2011')
                corrected_e  *= smearFactor
                corrected_et *= smearFactor
            else:
                ## Calibration correction in Data
                sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                scaleFactor = 1
                if corrected_e > 0:
                    scaleFactor = egammaER.applyEnergyCorrectionMeV(cl_eta, el.cl_phi, corrected_e, corrected_et, sys, 'ELECTRON') / corrected_e
                corrected_e  *= scaleFactor
                corrected_et *= scaleFactor


            ## Modify E and Et in transient D3PD
            el.cl_E  = corrected_e
            el.cl_et = corrected_et

            ## Modify the fourvector which is used for the rest of the tree filling
            el_eta = el.fourvect.Eta()
            el_phi = el.fourvect.Phi()
            newElfourvect = TLorentzVector()
            newElfourvect.SetPtEtaPhiE(corrected_et, el_eta, el_phi, corrected_e)
            setattr(el, 'fourvect', newElfourvect)
            

            ## Adjust MET accordingly
            if ((el.nSCTHits + el.nPixHits) < 4):
                px = raw_et * cos(el.cl_phi)
                py = raw_et * sin(el.cl_phi)
                corrected_px = corrected_et * cos(el.cl_phi)
                corrected_py = corrected_et * sin(el.cl_phi)
            else:
                px = raw_et * cos(el.trackphi)
                py = raw_et * sin(el.trackphi)
                corrected_px = corrected_et * cos(el.trackphi)
                corrected_py = corrected_et * sin(el.trackphi)

            event.MET_RefFinal_BDTMedium_etx += ( px - corrected_px ) * el.MET_BDTMedium_wpx[0]
            event.MET_RefFinal_BDTMedium_ety += ( py - corrected_py ) * el.MET_BDTMedium_wpy[0]
            event.MET_RefFinal_BDTMedium_sumet -= ( raw_et - corrected_et ) * el.MET_BDTMedium_wet[0]

        return True

    
from ROOT import CaloIsoCorrection
class ElectronIsoCorrection(EventFilter):
    """
    Correction for electron calorimeter isolation variables
    """
    
    def __init__(self, datatype, **kwargs):

        super(ElectronIsoCorrection, self).__init__(**kwargs)
        self.datatype = datatype

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
            

            if self.datatype == datasets.MC:
                newElectronEtCone20 = CaloIsoCorrection.GetPtNPVCorrectedIsolation(nPV,
                                                                                   E,
                                                                                   EtaS2,
                                                                                   EtaP,
                                                                                   cl_eta,
                                                                                   20,
                                                                                   True,
                                                                                   EtCone20)
            else:
                newElectronEtCone20 = CaloIsoCorrection.GetPtNPVCorrectedIsolation(nPV,
                                                                                   E,
                                                                                   EtaS2,
                                                                                   EtaP,
                                                                                   cl_eta,
                                                                                   20,
                                                                                   False,
                                                                                   EtCone20)
            el.Etcone20 = newEtCone20
        return True
        
    


#################################################
# Tau/Electron Misidentification correction
#################################################

def TauEfficiencySF(event, datatype):
    """
    Apply Tau Efficiency Scale Factor correction
    """

    # Apply only on MC
    if datatype != datasets.MC: return 1.0

    for tau in event.taus:
        #Correct 1p taus/ electron fake rate
        #https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012#Systematics_for_Electron_Misiden
        if tau.numTrack == 1 and tau.fourvect.Pt() > 20*GeV:
            nMC = event.mc_n
            for i in range(0, nMC):
                if abs(event.mc_pdgId[i]) == 11 and event.mc_pt[i] > 8*GeV:
                    if utils.dR(event.mc_eta[i], event.mc_phi[i], tau.eta, tau.phi) > 0.2: continue
                    if abs(tau.eta) < 1.37:
                        return 1.64
                    elif abs(tau.eta) < 1.52:
                        return 1.0
                    elif abs(tau.eta) < 2.0:
                        return 0.71
                    else:
                        return 2.9

    return 1.0


#################################################
# Muon Efficiency corrections
#################################################

## Scale factors for single muon triggers
from ROOT import Analysis
from ROOT import LeptonTriggerSF

leptonTriggerSF = LeptonTriggerSF(TrigMuonEfficiency.RESOURCE_PATH)

__cached__ = False

def getMuonSF(pileup):
    global __cached__
    if not __cached__:
        int_lum = pileup.getIntegratedLumiVector()
        __cached__ = Analysis.AnalysisMuonEfficiencyScaleFactors("STACO_CB", int_lum, "MeV", MuonEfficiencyCorrections.RESOURCE_PATH)
    return __cached__

def MuonSF(event, datatype, pileup_tool):
    """
    Apply Muon Efficiency correction for trigger and others
    """
    # Apply only on MC
    if datatype != datasets.MC: return 1.0

    # Weight
    w = 1.0

    #Get SF tool
    muonSF = getMuonSF(pileup_tool)

    for mu in event.muons:
        muon = TLorentzVector()
        muon.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, mu.m)
        w *= muonSF.scaleFactor(muon)

    return w

## Scale factor for lephad triggers
from ROOT import HSG4TriggerSF
def MuonLTTSF(muon, runNumber):

    sfTool = HSG4TriggerSF(HSG4.RESOURCE_PATH)
    return sfTool.getSFMuon(muon.fourvect, runNumber, 0)
        


#################################################
# Electron Efficiency corrections
#################################################

## Scale factors for single electron triggers
from ROOT import egammaSFclass

egammaSF = egammaSFclass()

def ElectronSF(event, datatype):

    ## Apply only on MC and embedded
    if datatype == datasets.DATA: return 1.0

    ## Weight
    weight = 1.0

    """
    Electron efficiency correction.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EfficiencyMeasurements
    """

    for iel, el in enumerate(event.electrons):
        if iel > 0:
            print 'WARNING in ElectronSF(): using more than 1 electron in efficiency correction'
        pt = el.fourvect.Pt()
        eta = el.cl_eta
        reconstructionEffCorr = (egammaSF.scaleFactor(eta, pt, 4, 0, 6, True)).first # track + reco eff SF
        identificationEffCorr = (egammaSF.scaleFactor(eta, pt, 7, 0, 6, True)).first # ID eff SF
        ## We don't use electrons in the forward region, so we don't need to apply SF for them.
        totalEffCorr = reconstructionEffCorr * identificationEffCorr
        weight *= totalEffCorr

    return weight


## Scale factor for single electron triggers
def LeptonSLTSF(event, datatype):
    """
    Trigger efficiency correction.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Electrons
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/TrigMuonEfficiency#LeptonTriggerSF_Tool
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
    
    trigSF = (leptonTriggerSF.GetTriggerSF(event.RunNumber, False, v_muTLVs, 1, v_elTLVs, 2)).first

    return trigSF

## Scale factor for lephad triggers
def ElectronLTTSF(electron, runNumber):

    sfTool = HSG4TriggerSF(HSG4.RESOURCE_PATH)
    return sfTool.getSFElec(electron.fourvect, runNumber, 0)



#################################################
# Special LTT corrections for embedded
#################################################
from ROOT import TauTriggerCorrections

def EmbedTauTriggerCorr(Tau, nvtx, runNumber):
    ttc = TauTriggerCorrections()
    weight = 1.0
    ttcPath = TTC.RESOURCE_PATH
    tauPt  = Tau.fourvect.Pt()
    tauEta = Tau.fourvect.Eta()
    
    if runNumber > 186755:
        status =  ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_wmcpara_EF_tau20_medium1.root'), '1P3P', 'BDTm')
        weight =  ttc.get3DMCEff(tauPt, tauEta , nvtx, 0)
        status =  ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_EF_tau20_medium1.root'))
        weight *= ttc.getSF(tauPt, 0)
    else:
        status = ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_EF_tau16_loose.root'))
        weight = ttc.getDataEff(tauPt, 0)

    return weight
            
        

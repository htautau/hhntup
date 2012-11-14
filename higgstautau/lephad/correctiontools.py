from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
import os

from externaltools.bundle_2012 import MuonIsolationCorrection
from externaltools import egammaAnalysisUtils
from externaltools import HSG4TriggerSF as HSG4

from ROOT import std, TLorentzVector, TFile
from higgstautau.mixins import *


#################################################
# Muon Pt Smearing
#################################################
__mmc_cached__  = False

def getMMC(year):
    global __mmc_cached__
    if not (__mmc_cached__):
        # if year == 2011:
        #     from externaltools.bundle_2011 import MuonMomentumCorrections as MMC2011
        #     __mmc_path_cached__ = MMC2011.RESOURCE_PATH
        # if year == 2012:
        #     from externaltools.bundle_2012 import MuonMomentumCorrections as MMC2012
        #     __mmc_path_cached__ = MMC2012.RESOURCE_PATH

        from externaltools import MuonMomentumCorrections as MMC
        from ROOT import MuonSmear
        if year == 2011:
            ## Prepare 2011 muon pt smearing tool
            __mmc_cached__ = MuonSmear.SmearingClass('Data11',
                                                     'staco',
                                                     'pT',
                                                     'Rel17',
                                                     MMC.RESOURCE_PATH)
            __mmc_cached__.UseScale(1)
            __mmc_cached__.UseImprovedCombine()
            __mmc_cached__.RestrictCurvatureCorrections(2.5)
            __mmc_cached__.FillScales('KC')

        if year == 2012:
            ## Prepare 2012 muon pt smearing tool
            __mmc_cached__ = MuonSmear.SmearingClass('Data12',
                                                     'staco',
                                                     'pT',
                                                     'Rel17.2_preliminary',
                                                     MMC.RESOURCE_PATH)
            __mmc_cached__.UseScale(1)
            __mmc_cached__.UseImprovedCombine()
            
    return __mmc_cached__

class MuonPtSmearing(EventFilter):
    """
    Smears the Muon Pt using the official ATLAS tool
    """

    def __init__(self, datatype, year, **kwargs):

        super(MuonPtSmearing, self).__init__(**kwargs)
        self.datatype = datatype

        self.year = year
        self.tool = getMMC(year)
        
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

            
            ## Seed with event number, reproducible smear for different analyses
            self.tool.SetSeed(event.EventNumber, i)

            ## Pass parameters, get smeared Pt
            pt_smear = -1

            if mu.isCombinedMuon:

                self.tool.Event(pt_ms, pt_id, pt, eta, charge)
                pt_smear = self.tool.pTCB()
                if self.year == 2012:
                    self.tool.Event(pt_standalone_forMET, eta, 'MS', charge)
                    pt_standalone_forMET = self.tool.pTMS()
                if self.year == 2011:
                    pt_standalone_forMET = pt_smear
                
            elif mu.isSegmentTaggedMuon:
                self.tool.Event(pt_id, eta, 'ID', charge)
                pt_smear = self.tool.pTID()
                if self.year == 2012:
                    self.tool.Event(pt_standalone_forMET, eta, 'MS', charge)
                    pt_standalone_forMET = self.tool.pTMS()
                if self.year == 2011:
                    pt_standalone_forMET = pt_smear

            else:
                self.tool.Event(pt_ms, eta, 'MS', charge)
                pt_smear = self.tool.pTMS()
                if self.year == 2012:
                    pt_standalone_forMET = pt_smear
                if self.year == 2011:
                    self.tool.Event(pt_standalone_forMET, eta, 'MS', charge)
                    pt_standalone_forMET = self.tool.pTMS()

            ## Adjust Pt in transient D3PD
            mu.pt = pt_smear

            ## Adjust Pt of the muon 4-vector
            newMufourvect = TLorentzVector()
            newMufourvect.SetPtEtaPhiM(pt_smear, eta, mu.fourvect.Phi(), mu.fourvect.M())
            setattr(mu, 'fourvect', newMufourvect)

            ## Adjust MET accordingly
            px = pt*cos(mu.phi)
            py = pt*sin(mu.phi)
            
            if self.year == 2011:
                px_smear = pt_standalone_forMET*cos(mu.phi)
                py_smear = pt_standalone_forMET*sin(mu.phi)

                event.MET_RefFinal_BDTMedium_etx += (px-px_smear)*mu.MET_BDTMedium_wpx[0]
                event.MET_RefFinal_BDTMedium_ety += (py-py_smear)*mu.MET_BDTMedium_wpy[0]
                event.MET_RefFinal_BDTMedium_sumet -= (pt-pt_smear)*mu.MET_BDTMedium_wet[0]

            if self.year == 2012:
                px_smear = pt_standalone_forMET*cos(mu.phi)
                py_smear = pt_standalone_forMET*sin(mu.phi)

                event.MET_RefFinal_STVF_etx += (px-px_smear)*mu.MET_STVF_wpx[0]
                event.MET_RefFinal_STVF_ety += (py-py_smear)*mu.MET_STVF_wpy[0]
                event.MET_RefFinal_STVF_sumet -= (pt-pt_smear)*mu.MET_STVF_wet[0]

        return True


#################################################
# Muon isolation correction
#################################################
from ROOT import CorrectCaloIso

class MuonIsoCorrection(EventFilter):
    """
    Corrects the muon isolation in 2012
    """

    def __init__(self, datatype, year, **kwargs):

        super(MuonIsoCorrection, self).__init__(**kwargs)
        self.datatype = datatype
        self.year = year
        self.tool = CorrectCaloIso()

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
# Electron Energy Rescaling
#################################################
from ROOT import eg2011

class EgammaERescaling(EventFilter):
    """
    Rescales/smeares the electron energy using the official ATLAS tool.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EnergyScaleResolutionRecommendations
    """

    def __init__(self, datatype, year, **kwargs):

        super(EgammaERescaling, self).__init__(**kwargs)
        self.datatype = datatype
        self.year = year
        self.tool = eg2011.EnergyRescaler()
        if self.year == 2011:
            self.tool.useDefaultCalibConstants('2011')
        if self.year == 2012:
            self.tool.useDefaultCalibConstants('2012')

    def passes(self, event):

        for el in event.electrons:

            ## Seed with event number, reproducible smear for different analyses
            self.tool.SetRandomSeed(int(5*event.RunNumber+event.EventNumber+(getattr(el,'fourvect').Phi()+pi)*1000000.))

            raw_e  = el.cl_E
            raw_et = getattr(el,'fourvect').Pt()
            corrected_e  = raw_e
            corrected_et = raw_et
            cl_eta = el.cl_eta
            trk_eta = el.tracketa
            cl_phi = el.cl_phi

            if self.year == 2011:
                ## Calibration for electrons in the transition region in Data and MC
                scaleFactorForTransitionRegion = self.tool.applyMCCalibrationMeV(cl_eta, raw_et, 'ELECTRON')
                corrected_e  *= scaleFactorForTransitionRegion
                corrected_et *= scaleFactorForTransitionRegion

                if self.datatype == datasets.MC or self.datatype == datasets.EMBED:
                    ## Smearing correction in MC
                    sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                    smearFactor = self.tool.getSmearingCorrectionMeV(cl_eta, corrected_e, sys, False, '2011')
                    corrected_e  *= smearFactor
                    corrected_et *= smearFactor
                else:
                    ## Calibration correction in Data
                    sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                    scaleFactor = 1
                    if corrected_e > 0:
                        scaleFactor = self.tool.applyEnergyCorrectionMeV(cl_eta, cl_phi, corrected_e, corrected_et, sys, 'ELECTRON') / corrected_e
                    corrected_e  *= scaleFactor
                    corrected_et *= scaleFactor

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

            if self.year == 2011:
                event.MET_RefFinal_BDTMedium_etx += ( px - corrected_px ) * el.MET_BDTMedium_wpx[0]
                event.MET_RefFinal_BDTMedium_ety += ( py - corrected_py ) * el.MET_BDTMedium_wpy[0]
                event.MET_RefFinal_BDTMedium_sumet -= ( raw_et - corrected_et ) * el.MET_BDTMedium_wet[0]
            if self.year == 2012:
                event.MET_RefFinal_STVF_etx += ( px - corrected_px ) * el.MET_STVF_wpx[0]
                event.MET_RefFinal_STVF_ety += ( py - corrected_py ) * el.MET_STVF_wpy[0]
                event.MET_RefFinal_STVF_sumet -= ( raw_et - corrected_et ) * el.MET_STVF_wet[0]
                

        return True


#################################################
# Electron Isolation correction
#################################################
    
from ROOT import CaloIsoCorrection
class ElectronIsoCorrection(EventFilter):
    """
    Correction for electron calorimeter isolation variables
    """
    
    def __init__(self, datatype, year, **kwargs):

        self.datatype = datatype
        self.year = year        
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

            if self.year == 2012:
                newTopoEtCone20 = CaloIsoCorrection.GetPtEDCorrectedTopoIsolation(el.ED_median,
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

from rootpy.io import open as ropen
HERE = os.path.dirname(os.path.abspath(__file__))
MuonIsoCorrFile = ropen(os.path.join(HERE, 'SF_2D_UptoE5.root'))
Corr2D = MuonIsoCorrFile.Get('SF_2D_PtVsNvx')

def MuonIsoEffCorrection(event, datatype, year):
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

    b = Corr2D.FindBin(mu.fourvect.Pt()/1000, hsg3_nvtx)
    weight = Corr2D.GetBinContent(b)

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
from externaltools import MuonEfficiencyCorrections
from ROOT import Analysis

__cached__ = False

def getMuonSF(pileup, year):
    global __cached__
    if not __cached__:
        if year == 2011:
            int_lum = pileup.getIntegratedLumiVector()
            __cached__ = Analysis.AnalysisMuonEfficiencyScaleFactors('STACO_CB',
                                                                     int_lum,
                                                                     'MeV',
                                                                     MuonEfficiencyCorrections.RESOURCE_PATH)
        if year == 2012:
            __cached__ = Analysis.AnalysisMuonConfigurableScaleFactors('',
                                                                       MuonEfficiencyCorrections.RESOURCE_PATH + '/STACO_CB_2012_SF.txt',
                                                                       'MeV',
                                                                       Analysis.AnalysisMuonConfigurableScaleFactors.AverageOverPeriods)
            __cached__.setRunInterval(200804, 210308)
            __cached__.Initialise()
            
    return __cached__

def MuonSF(event, datatype, pileup_tool, year):
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
    muonSF = getMuonSF(pileup_tool, year)

    if year == 2011:
        sf = muonSF.scaleFactor(muon)
        sf_err_stat = muonSF.scaleFactorUncertainty(muon)
        sf_err_syst = muonSF.scaleFactorSystematicUncertainty(muon)
        errup = sqrt(sf_err_stat**2 + sf_err_syst**2)
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
from ROOT import HSG4TriggerSF
sfTool = HSG4TriggerSF(HSG4.RESOURCE_PATH)

def MuonLTTSF(event, datatype, year, pileup_tool):
    """
    Returns nominal, +1sigma, -1sigma
    The -1sigma and +1sigma weights are designed to be applied on top of everything else
    """

    sf, errup, errdown = 1.0, 0.0, 0.0

    mu = event.muons[0].fourvect
    runNumber = event.RunNumber

    if year == 2011:
        if datatype == datasets.MC:
            sf      = sfTool.getSFMuon(mu, runNumber, 0)
            errup   = sfTool.getSFMuon(mu, runNumber, 1)
            errdown = sfTool.getSFMuon(mu, runNumber, -1)
        else:
            sf      = sfTool.getDataEffMuon(mu, runNumber, 0)
            errup   = sfTool.getDataEffMuon(mu, runNumber, 1)
            errdown = sfTool.getDataEffMuon(mu, runNumber, -1)


    elif year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(runNumber)
        
        if datatype == datasets.MC:
            sf      = sfTool.getSFMuon(mu, randomRunNumber, 0)
            errup   = sfTool.getSFMuon(mu, randomRunNumber, 1)
            errdown = sfTool.getSFMuon(mu, randomRunNumber, -1)
        else:
            sf      = sfTool.getDataEffMuon(mu, randomRunNumber, 0)
            errup   = sfTool.getDataEffMuon(mu, randomRunNumber, 1)
            errdown = sfTool.getDataEffMuon(mu, randomRunNumber, -1)

    if sf > 0.0:
        errup   = errup / sf
        errdown = errdown / sf

    return sf, errup, errdown
        
        


#################################################
# Electron Efficiency corrections
#################################################

## Scale factors for single electron triggers
from ROOT import egammaSFclass
egammaSF = egammaSFclass()

def ElectronSF(event, datatype, pileup_tool, year):

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
        reconstructionEffCorr = egammaSF.scaleFactor(eta, pt, 4, 0, 6, True) # track + reco eff SF
        identificationEffCorr = egammaSF.scaleFactor(eta, pt, 7, 0, 6, True) # ID eff SF
    if year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(event.RunNumber)
        reconstructionEffCorr = egammaSF.scaleFactor(eta, pt, 4, 0, 8, 1, randomRunNumber) # track + reco eff SF
        identificationEffCorr = egammaSF.scaleFactor(eta, pt, 7, 0, 8, 1, randomRunNumber) # ID eff SF
    ## We don't use electrons in the forward region, so we don't need to apply SF for them.
    sf = reconstructionEffCorr.first * identificationEffCorr.first
    reconstructionError = reconstructionEffCorr.second
    identificationError = identificationEffCorr.second
    errup = sqrt(reconstructionError**2 + identificationError**2)
    errdown = errup

    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown




#################################################
# Lepton trigger scale factors
#################################################
__ltsf_cached__ = False

def getLTSF(year):
    global __ltsf_cached__
    if not __ltsf_cached__:
        if year == 2011:
            from externaltools.bundle_2011 import TrigMuonEfficiency
            from ROOT import LeptonTriggerSF
            __lstf_cached__ = LeptonTriggerSF(TrigMuonEfficiency.RESOURCE_PATH)
        if year == 2012:
            from externaltools.bundle_2012 import TrigMuonEfficiency
            from ROOT import LeptonTriggerSF
            __ltsf_cached__ = LeptonTriggerSF(2012, TrigMuonEfficiency.RESOURCE_PATH, 'muon_trigger_sf_2012_AtoE.root')
            
    return __ltsf_cached__


## Scale factor for single lepton triggers
def LeptonSLTSF(event, datatype, pileup_tool, year):

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
    if datatype == datasets.DATA: sf, errup, errdown

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

    randomRunNumber = pileup_tool.GetRandomRunNumber(event.RunNumber)

    leptonTriggerSF = getLTSF(year)

    trigSF = leptonTriggerSF.GetTriggerSF(randomRunNumber, False, v_muTLVs, 1, v_elTLVs, 2, 0)
    sf = trigSF.first
    errup = trigSF.second
    errdown = errup

    if sf > 0.0:
        errup   = (sf + errup) / sf
        errdown = (sf - errdown) / sf

    return sf, errup, errdown



## Scale factor for lephad triggers
def ElectronLTTSF(event, datatype, year, pileup_tool):
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
            sf      = sfTool.getSFElec(e, runNumber, 0)
            errup   = sfTool.getSFElec(e, runNumber, 1)
            errdown = sfTool.getSFElec(e, runNumber, -1)
        else:
            sf      = sfTool.getDataEffElec(e, runNumber, 0)
            errup   = sfTool.getDataEffElec(e, runNumber, 1)
            errdown = sfTool.getDataEffElec(e, runNumber, -1)

    elif year == 2012:
        randomRunNumber = pileup_tool.GetRandomRunNumber(runNumber)
        
        if datatype == datasets.MC:
            sf      = sfTool.getSFElec(e, randomRunNumber, 0)
            errup   = sfTool.getSFElec(e, randomRunNumber, 1)
            errdown = sfTool.getSFElec(e, randomRunNumber, -1)
        else:
            sf      = sfTool.getDataEffElec(e, randomRunNumber, 0)
            errup   = sfTool.getDataEffElec(e, randomRunNumber, 1)
            errdown = sfTool.getDataEffElec(e, randomRunNumber, -1)

    if sf > 0.0:
        errup   = errup / sf
        errdown = errdown / sf

    return sf, errup, errdown


    

#################################################
# Special LTT corrections for tau trigger
#################################################

__ttc_cached__  = False
__path_cached__ = False

def getTTC(year):
    global __ttc_cached__
    global __path_cached__
    if not (__ttc_cached__ or __path_cached__):
        if year == 2011:
            from externaltools.bundle_2011 import TauTriggerCorrections as TTC2011
            __path_cached__ = TTC2011.RESOURCE_PATH
        if year == 2012:
            from externaltools.bundle_2012 import TauTriggerCorrections as TTC2012
            __path_cached__ = TTC2012.RESOURCE_PATH

    from ROOT import TauTriggerCorrections
    __ttc_cached__ = TauTriggerCorrections()
    return __ttc_cached__, __path_cached__


def TauTriggerSF(Tau, nvtx, runNumber, year, lep, pileup_tool):
    
    sf, errup, errdown = 1.0, 0.0, 0.0
    ttc, ttcPath = getTTC(year)
    tauPt  = Tau.fourvect.Pt()
    tauEta = Tau.fourvect.Eta()

    if year == 2011:
        if runNumber > 186755:
            status    =  ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_wmcpara_EF_tau20_medium1.root'), '1P3P', 'BDTm')
            sf    =  ttc.get3DMCEff(tauPt, tauEta , nvtx, 0)
            errup    += ttc.get3DMCEff(tauPt, tauEta , nvtx, 1)**2
            errdown  += ttc.get3DMCEff(tauPt, tauEta , nvtx, -1)**2
            status    =  ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_EF_tau20_medium1.root'))
            sf   *= ttc.getSF(tauPt, 0)
            errup    += ttc.getSF(tauPt, 1)**2
            errdown  += ttc.getSF(tauPt, -1)**2
            errup = sqrt(errup)
            errdown = sqrt(errdown)
        else:
            status  = ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_EF_tau16_loose.root'))
            sf  = ttc.getDataEff(tauPt, 0)
            errup   = ttc.getDataEff(tauPt, 1)
            errdown = ttc.getDataEff(tauPt, -1)

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
            ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_EF_tau20_medium1.root'))
            sf  = ttc.getSF(tauPt, tauEta, 0, period, nprong, 'BDTm', 'EVm')
            errup   = ttc.getSF(tauPt, tauEta, 1, period, nprong, 'BDTm', 'EVm')
            errdown = ttc.getSF(tauPt, tauEta, -1, period, nprong, 'BDTm', 'EVm')

        if lep == 'e':
            ttc.loadInputFile(os.path.join(ttcPath, 'triggerSF_EF_tau20Ti_medium1.root'))
            sf  = ttc.getSF(tauPt, tauEta, 0, period, nprong, 'BDTm', 'EVm')
            errup   = ttc.getSF(tauPt, tauEta, 1, period, nprong, 'BDTm', 'EVm')
            errdown = ttc.getSF(tauPt, tauEta, -1, period, nprong, 'BDTm', 'EVm')

    if sf > 0.0:
        errup   = errup / sf
        errdown = errdown / sf

    return sf, errup, errdown
            
        

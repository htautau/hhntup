from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
from external.Muons     import muonSmear, getMuonSF, muonTriggerSF
from external.Electrons import egammaER
from external.ggF       import getggFTool
from ROOT import *
from higgstautau.mixins import *


#################################################
# Muon Pt Smearing
#################################################

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

        if self.datatype != datasets.MC: return True

        for mu in event.muons:

            #Obtain parameters for correction
            charge = mu.charge
            eta    = mu.eta
            pt     = mu.pt
            pt_ms  = sin(mu.ms_theta)/abs(mu.ms_qoverp)
            pt_id  = sin(mu.id_theta)/abs(mu.id_qoverp)

            #Seed with event number, reproducible smear for different analyses
            muonSmear.SetSeed(event.EventNumber)

            #Pass parameters, get smeared Pt
            muonSmear.Event(pt_ms, pt_id, pt, eta, charge)
            pt_smear = -1

            if mu.isCombinedMuon:
                pt_smear = muonSmear.pTCB()
            else:
                pt_smear = muonSmear.pTID()

            #Adjust Pt in transient D3PD
            mu.pt = pt_smear

            #Adjust MET accordingly
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

            # Seed with event number, reproducible smear for different analyses
            egammaER.SetRandomSeed(int(5*event.RunNumber+event.EventNumber+(getattr(el,'fourvect').Phi()+pi)*1000000.))

            raw_e  = el.cl_E
            raw_et = getattr(el,'fourvect').Pt()
            corrected_e  = raw_e
            corrected_et = raw_et

            # Calibration for electrons in the transition region in Data and MC
            scaleFactorForTransitionRegion = egammaER.applyMCCalibrationMeV(el.cl_eta, raw_et, 'ELECTRON')
            corrected_e  *= scaleFactorForTransitionRegion
            corrected_et *= scaleFactorForTransitionRegion

            if self.datatype == datasets.MC:
                # Smearing correction in MC
                sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                smearFactor = egammaER.getSmearingCorrectionMeV(el.cl_eta, corrected_e, sys, False, '2011')
                corrected_e  *= smearFactor
                corrected_et *= smearFactor
            else:
                # Calibration correction in Data
                sys = 0 # 0: nominal, 1: -1sigma, 2: +1sigma
                scaleFactor = 1
                if corrected_e > 0:
                    scaleFactor = egammaER.applyEnergyCorrectionMeV(el.cl_eta, el.cl_phi, corrected_e, corrected_et, sys, 'ELECTRON') / corrected_e
                corrected_e  *= scaleFactor
                corrected_et *= scaleFactor

 
            # Modify E and Et in transient D3PD
            el.cl_E  = corrected_e
            el.cl_et = corrected_et
            
            # Adjust MET accordingly
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
                    if tau.eta < 1.37:
                        return 1.64
                    elif tau.eta < 1.52:
                        return 1.0
                    elif tau.eta < 2.0:
                        return 0.71
                    else:
                        return 2.9

    return 1.0


#################################################
# Muon Efficiency corrections
#################################################

def MuonSF(event, datatype, pileup_tool):
    """
    Apply Muon Efficiency correction for trigger and others
    """
    # Apply only on MC
    if datatype != datasets.MC: return 1.0

    # Weight
    w = 1.0

    # Load muonSF tool
    muonSF = getMuonSF(pileup_tool)

    #Store electrons and muons for the trigger efficiency tool
    std_muons     = std.vector(TLorentzVector)()
    std_electrons = std.vector(TLorentzVector)()

    for mu in event.muons:
        muon = TLorentzVector()
        muon.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, mu.m)
        std_muons.push_back(muon)
        w *= muonSF.scaleFactor(muon)

    for e in event.electrons:
        electron = TLorentzVector()
        cl_eta = e.cl_eta
        trk_eta = e.tracketa
        cl_Et = e.cl_E/cosh(trk_eta)
        electron.SetPtEtaPhiM(cl_Et, cl_eta, e.cl_phi, e.m)
        std_electrons.push_back(electron)

    trigSF = muonTriggerSF.GetTriggerSF(event.RunNumber, false, std_muons, 1, std_electrons, 2)
    w *= trigSF.first

    return w
    

#################################################
# Electron Efficiency corrections
#################################################

def ElectronSF(event, datatype, pileupTool):

    ## Apply only on MC
    if datatype != datasets.MC: return 1.0

    ## Weight
    weight = 1.0

    """
    Electron efficiency correction.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/EfficiencyMeasurements
    """

    ## Load muonSF tool
    muonSF = getMuonSF(pileupTool)

    ## Store electrons and muons for the trigger efficiency tool
    v_muTLVs = std.vector(TLorentzVector)()
    v_elTLVs = std.vector(TLorentzVector)()

    for mu in event.muons:
        print "There should not be a muon"
        muTLV = TLorentzVector()
        muTLV.SetPtEtaPhiM(mu.pt, mu.eta, mu.phi, mu.m)
        v_muTLVs.push_back(muTLV)
        weight *= muonSF.scaleFactor(muTLV)

    for iel, el in enumerate(event.electrons):
        if iel > 0:
            print 'WARNING in ElectronSF(): using more than 1 electron in efficiency correction'
        elTLV = TLorentzVector()
        pt = getattr(el,'fourvect').Pt()
        eta = el.cl_eta
        phi = el.cl_phi
        E = el.cl_E
        elTLV.SetPtEtaPhiE(pt,eta,phi,E)
        v_elTLVs.push_back(elTLV)
        reconstructionEffCorr = (egammaSF.scaleFactor(el.cl_eta, getattr(el,'fourvect').Pt(), 4, 0, 6, True)).first # track + reco eff SF
        identificationEffCorr = (egammaSF.scaleFactor(el.cl_eta, getattr(el,'fourvect').Pt(), 7, 0, 6, True)).first # ID eff SF
        ## We don't use electrons in the forward region, so we don't need to apply SF for them.
        totalEffCorr = reconstructionEffCorr * identificationEffCorr
        weight *= totalEffCorr

    """
    Trigger efficiency correction.
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Summer#Electrons
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/TrigMuonEfficiency#LeptonTriggerSF_Tool
    """
    if v_muTLVs.size() > 0: print 'WARNING in ElectronSF(): vector of muons is NOT empty'
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
    pileupTool.SetRandomSeed(314159 + event.mc_channel_number*2718 + event.EventNumber)
    randomRunNumber = pileupTool.GetRandomRunNumber(event.RunNumber)
    trigSF = (leptonTriggerSF.GetTriggerSF(randomRunNumber, False, v_muTLVs, 1, v_elTLVs, 2)).first
    weight *= trigSF

    return weight


#################################################
#ggF reweighting
#################################################

def ggFreweighting(event, dataname):
    """
    Reweight the ggF samples
    """

    if dataname.find('PowHegPythia_ggH') == -1: return 1.0

    #Extract mass value
    mass = int(dataname.lstrip('PowHegPythia_ggH').rstrip('_tautaulh.mc11c'))

    #Get corresponding ggF tool setting
    ggFTool = getggFTool(mass)

    #Find the Higgs particle in mc
    nMC = event.mc_n
    pt = 0
    for i in range(0, nMC):
        if event.mc_pdgId[i] == 25 and event.mc_status[i] != 3:
            pt = event.mc_pt[i]/1000

    if pt > 0:
        return ggFTool.getWeight(pt)
    else:
        return 1.0

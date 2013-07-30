"""
This module defines the output branches in the final ntuple
as TreeModels.
"""


from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy.tree.treetypes import *

from atlastools.utils import et2pt
from atlastools.units import GeV
from atlastools import utils

from ..models import MatchedObject, TrueTau

import math
import ROOT

##################################################
class RecoTau(TreeModel):

    BDTJetScore = FloatCol()
    BDTEleScore = FloatCol()

    JetBDTSigLoose = BoolCol()
    JetBDTSigMedium = BoolCol()
    JetBDTSigTight = BoolCol()
    EleBDTtight = BoolCol()
    isTrueLep = BoolCol()

    numTrack = IntCol()
    charge = IntCol()

    fourvect = LorentzVector


    
##################################################
class RecoLepton(TreeModel):

    isolated = BoolCol()
    charge = IntCol()
    fourvect = LorentzVector
    leptype = BoolCol()


##################################################
class RecoMET(TreeModel):

    MET_vect = Vector2
    MET      = FloatCol()

    
##################################################
class EventVariables(TreeModel):

    dilep_veto = BoolCol()
    dilep_control = BoolCol()
    LTT = BoolCol()
    SLT = BoolCol()

    FF = FloatCol()
    FF_up = FloatCol()
    FF_down = FloatCol()

    is_tau = BoolCol()

    subsample1 = BoolCol()
    subsample2 = BoolCol()

    subsample2_1 = BoolCol()
    subsample2_2 = BoolCol()

    category_vbf_train = BoolCol()
    category_vbf_test = BoolCol()
    category_boosted = BoolCol()
    category_1j = BoolCol()
    category_0j = BoolCol()
    
    tau_x = FloatCol()
    lep_x = FloatCol()
    tau_x_lep_x = FloatCol()
    mass_mmc_tau_lep = FloatCol()
    pt_mmc_tau_lep = FloatCol()
    met_mmc_tau_lep = FloatCol()
    mass_vis_tau_lep = FloatCol()
    pt_ratio_tau_lep = FloatCol()
    dphi_met_lep = FloatCol()
    pt_vector_sum_all = FloatCol()

    numJets = IntCol()
    numJets30 = IntCol()
    numJets50 = IntCol()
    btag = BoolCol()

    leadjet_fourvect = LorentzVector
    subleadjet_fourvect = LorentzVector
    
    sumPt = FloatCol()
    cutflow = IntCol()

    charge_product_tau_lep = IntCol()
    mass_transverse_met_lep = FloatCol()
    mass_transverse_met_tau = FloatCol()
    ddr_tau_lep = FloatCol()
    dr_tau_lep = FloatCol()
    resonance_pt_tau_lep = FloatCol()

    mass_j1_j2 = FloatCol()
    eta_product_j1_j2 = FloatCol()
    eta_delta_j1_j2 = FloatCol()
    tau_centrality_j1_j2 = FloatCol()
    lep_centrality_j1_j2 = FloatCol()
    tau_lep_centrality_j1_j2 = FloatCol()
    met_phi_centrality = FloatCol()
    leadJetPt = FloatCol()
    sphericity = FloatCol()

    weight = FloatCol()
    nvtx = IntCol()

    true_higgs_mass = FloatCol()
    true_dphi_resonance_dijet = FloatCol()


##################################################
class SysWeights(TreeModel):

    sys_tau_ESF_UP = FloatCol()
    sys_tau_ESF_DOWN = FloatCol()
    sys_tau_IDSF_UP = FloatCol()
    sys_tau_IDSF_DOWN = FloatCol()
    sys_tau_TRIGSF_UP = FloatCol()
    sys_tau_TRIGSF_DOWN = FloatCol()

    sys_mu_TRIGSF_UP = FloatCol()
    sys_mu_TRIGSF_DOWN = FloatCol()
    sys_mu_ISOSF_UP = FloatCol()
    sys_mu_ISOSF_DOWN = FloatCol()
    sys_mu_EFFSF_UP = FloatCol()
    sys_mu_EFFSF_DOWN = FloatCol()

    sys_e_TRIGSF_UP = FloatCol()
    sys_e_TRIGSF_DOWN = FloatCol()
    sys_e_EFFSF_UP = FloatCol()
    sys_e_EFFSF_DOWN = FloatCol()

    w_mc           = FloatCol()
    w_pileup       = FloatCol()
    w_lumi         = FloatCol()
    w_tauesf       = FloatCol()
    w_tauidsf      = FloatCol()
    w_leptonsf     = FloatCol()
    w_leptontrigsf = FloatCol()
    w_tautriggersf = FloatCol()
    w_muonisosf    = FloatCol()
    w_ggf          = FloatCol()


##################################################
class RecoTauLepBlock((RecoTau).prefix('tau_') + (RecoLepton).prefix('lep_')):

    @classmethod
    def set(cls, event, tree, tau, lep, leptype, isMC=True, year=2011):
        """
        Misc variables
        """
        tree.mass_vis_tau_lep = (tau.fourvect + lep.fourvect).M()
        tree.charge_product_tau_lep = tau.charge * lep.charge
        tree.pt_ratio_tau_lep = lep.fourvect.Pt()/tau.fourvect.Pt()
        tree.pt_balance_tau_lep = (lep.fourvect.Pt() - tau.fourvect.Pt())/(tau.fourvect.Pt() + lep.fourvect.Pt())

        #Set tau variables
        setattr(tree, 'tau_BDTJetScore', tau.BDTJetScore)
        setattr(tree, 'tau_BDTEleScore', tau.BDTEleScore)
        setattr(tree, 'tau_JetBDTSigLoose', tau.JetBDTSigLoose)
        setattr(tree, 'tau_JetBDTSigMedium', tau.JetBDTSigMedium)
        setattr(tree, 'tau_JetBDTSigTight', tau.JetBDTSigTight)
        setattr(tree, 'tau_EleBDTtight', tau.EleBDTTight)
        setattr(tree, 'tau_numTrack', tau.numTrack)
        setattr(tree, 'tau_charge', tau.charge)
        getattr(tree, 'tau_fourvect').set_from(tau.fourvect)

        #Set lepton variables
        setattr(tree, 'lep_charge', lep.charge)
        getattr(tree, 'lep_fourvect').set_from(lep.fourvect)
        tree.lep_leptype = leptype

        #Calculate muon isolation
        if leptype == 0: # Is a muon
            muon_pt = lep.fourvect.Pt()
            track_iso = lep.ptcone40/muon_pt <= 0.06
            calo_iso = False
            if year == 2011:
                calo_iso  = lep.etcone20/muon_pt <= 0.04
            if year == 2012:
                calo_iso  = lep.etcone20/muon_pt <= 0.06
            setattr(tree, 'lep_isolated', (track_iso and calo_iso))

        #Calculate electron isolation
        if leptype == 1: # Is an electron
            track_iso = ( lep.ptcone40 / lep.cl_et < 0.06 )
            calo_iso = False
            if year == 2011:
                calo_iso  = ( lep.Etcone20 / lep.cl_et < 0.08 )
            if year == 2012:
                calo_iso  = lep.topoEtcone20 / lep.cl_et < 0.06
            setattr(tree, 'lep_isolated', (track_iso and calo_iso))
            

        #Find out if the tau is matched to an electron
        TauIsLep = False
        
        if tau.numTrack == 1 and tau.fourvect.Pt() > 20*GeV and isMC:
            nMC = event.mc_n
            for i in range(0, nMC):
                if abs(event.mc_pdgId[i]) == 11 and event.mc_pt[i] > 8*GeV:
                    if utils.dR(event.mc_eta[i], event.mc_phi[i], tau.eta, tau.phi) < 0.2: TauIsLep = True

        setattr(tree, "tau_isTruthLep", TauIsLep)

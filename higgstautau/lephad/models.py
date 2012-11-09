"""
This module defines the output branches in the final ntuple
as TreeModels.
"""


from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy.types import *

from atlastools.utils import et2pt
from atlastools.units import GeV
from atlastools import utils

from ..models import MatchedObject, TrueTau

import math
import ROOT


class RecoTau(TreeModel):

    BDTJetScore = FloatCol()
    BDTEleScore = FloatCol()

    JetBDTSigLoose = BoolCol()
    JetBDTSigMedium = BoolCol()
    JetBDTSigTight = BoolCol()
    EleBDTtight = BoolCol()
    isTruthEl = BoolCol()

    numTrack = IntCol()
    charge = IntCol()

    fourvect = LorentzVector


class RecoLepton(TreeModel):

    isolated = BoolCol()
    charge = IntCol()
    fourvect = LorentzVector
    leptype = BoolCol()

    BDTJetLoose  = IntCol()
    BDTJetMedium = IntCol()
    BDTJetTight  = IntCol()

    BDTEleLoose  = IntCol()
    BDTEleMedium = IntCol()
    BDTEleTight  = IntCol()


class RecoMET(TreeModel):

    MET_vect = Vector2
    MET_sig  = FloatCol()
    MET      = FloatCol()


class RecoJet(TreeModel):

    fourvect = LorentzVector
    
class EventVariables(TreeModel):

    mass_collinear_tau_lep = FloatCol()
    tau_x = FloatCol()
    lep_x = FloatCol()
    mass_mmc_tau_lep = FloatCol()
    pt_mmc_tau_lep = FloatCol()
    met_mmc_tau_lep = FloatCol()
    mass_vis_tau_lep = FloatCol()
    mass2_vis_tau_lep = FloatCol()
    cos_theta_tau_lep = FloatCol()
    pt_ratio_tau_lep = FloatCol()
    dphi_met_lep = FloatCol()

    numJets = IntCol()
    numJets30 = IntCol()
    numJets50 = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    truthjet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    jet_btag  = ROOT.vector('float')
    numVertices = IntCol()
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
    met_phi_centrality = FloatCol()
    tau_j1_j2_phi_centrality = FloatCol()

    leadJetPt = FloatCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    weight = FloatCol()

    nvtx = IntCol()


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
    sys_mu_SLTSF_UP = FloatCol()
    sys_mu_SLTSF_DOWN = FloatCol()
    sys_mu_LTTSF_UP = FloatCol()
    sys_mu_LTTSF_DOWN = FloatCol()

    sys_e_TRIGSF_UP = FloatCol()
    sys_e_TRIGSF_DOWN = FloatCol()
    sys_e_EFFSF_UP = FloatCol()
    sys_e_EFFSF_DOWN = FloatCol()
    sys_e_SLTSF_UP = FloatCol()
    sys_e_SLTSF_DOWN = FloatCol()
    sys_e_LTTSF_UP = FloatCol()
    sys_e_LTTSF_DOWN = FloatCol()


class RecoTauLepBlock((RecoTau).prefix('tau_') + (RecoLepton).prefix('lep_')):

    @classmethod
    def set(cls, event, tree, tau, lep, leptype, isMC=True, year=2011):
        """
        Misc variables
        """
        tree.mass_vis_tau_lep = utils.Mvis(tau.Et, tau.phi, lep.fourvect.Pt(), lep.fourvect.Phi())
        tree.mass2_vis_tau_lep = (tau.fourvect + lep.fourvect).M()
        theta_tau_lep = tau.fourvect.Vect().Angle(lep.fourvect.Vect())
        tree.cos_theta_tau_lep = math.cos(theta_tau_lep)
        tree.charge_product_tau_lep = tau.charge * lep.charge
        tree.pt_ratio_tau_lep = lep.fourvect.Pt()/tau.fourvect.Pt()

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
        TauIsEl = False
        
        if tau.numTrack == 1 and tau.fourvect.Pt() > 20*GeV and isMC:
            nMC = event.mc_n
            for i in range(0, nMC):
                if abs(event.mc_pdgId[i]) == 11 and event.mc_pt[i] > 8*GeV:
                    if utils.dR(event.mc_eta[i], event.mc_phi[i], tau.eta, tau.phi) < 0.2: TauIsEl = True

        setattr(tree, "tau_isTruthEl", TauIsEl)

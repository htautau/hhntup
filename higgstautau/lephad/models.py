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


class RecoMuon(TreeModel):

    isolated = BoolCol()
    charge = IntCol()
    fourvect = LorentzVector


class RecoElectron(TreeModel):

    isolated = BoolCol()
    charge = IntCol()
    fourvect = LorentzVector

class RecoLepton(TreeModel):

    isolated = BoolCol()
    charge = IntCol()
    fourvect = LorentzVector
    leptype = BoolCol()


class RecoMET(TreeModel):

    MET_vect = Vector2
    MET_sig  = FloatCol()
    MET      = FloatCol()


class RecoJet(TreeModel):

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector

    jvtxf = FloatCol()


class TauElectronEventVariables(TreeModel):

    mass_collinear_tau_electron = FloatCol()
    tau_x = FloatCol()
    electron_x = FloatCol()
    mass_mmc_tau_electron = FloatCol()
    pt_mmc_tau_electron = FloatCol()
    met_mmc_tau_electron = FloatCol()
    mass_vis_tau_electron = FloatCol()
    mass2_vis_tau_electron = FloatCol()
    theta_tau_electron = FloatCol()
    cos_theta_tau_electron = FloatCol()
    pt_ratio_tau_electron = FloatCol()
    dphi_met_electron = FloatCol()

    numJets = IntCol()
    numJets30 = IntCol()
    numJets50 = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    jet_btag  = ROOT.vector('float')
    numVertices = IntCol()
    HT = FloatCol()
    cutflow = IntCol()

    charge_product_tau_electron = IntCol()
    mass_transverse_met_electron = FloatCol()
    mass_transverse_met_tau = FloatCol()
    ddr_tau_electron = FloatCol()
    dr_tau_electron = FloatCol()
    higgs_pt = FloatCol()

    mass_j1_j2 = FloatCol()
    eta_product_j1_j2 = FloatCol()
    eta_delta_j1_j2 = FloatCol()
    tau_centrality_j1_j2 = FloatCol()
    electron_centrality_j1_j2 = FloatCol()
    met_phi_centrality = FloatCol()
    tau_j1_j2_phi_centrality = FloatCol()
    
    neff_pt = FloatCol()
    mass_all_jets = FloatCol()

    leadJetPt = FloatCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    weight = FloatCol()

    nvtx = IntCol()


class TauMuonEventVariables(TreeModel):

    mass_collinear_tau_muon = FloatCol()
    tau_x = FloatCol()
    muon_x = FloatCol()
    mass_mmc_tau_muon = FloatCol()
    pt_mmc_tau_muon = FloatCol()
    met_mmc_tau_muon = FloatCol()
    mass_vis_tau_muon = FloatCol()
    mass2_vis_tau_muon = FloatCol()
    theta_tau_muon = FloatCol()
    cos_theta_tau_muon = FloatCol()
    pt_ratio_tau_muon = FloatCol()
    dphi_met_muon = FloatCol()

    numJets = IntCol()
    numJets30 = IntCol()
    numJets50 = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    jet_btag  = ROOT.vector('float')
    numVertices = IntCol()
    HT = FloatCol()
    cutflow = IntCol()

    charge_product_tau_muon = IntCol()
    mass_transverse_met_muon = FloatCol()
    mass_transverse_met_tau = FloatCol()
    ddr_tau_muon = FloatCol()
    dr_tau_muon = FloatCol()
    higgs_pt = FloatCol()

    mass_j1_j2 = FloatCol()
    eta_product_j1_j2 = FloatCol()
    eta_delta_j1_j2 = FloatCol()
    tau_centrality_j1_j2 = FloatCol()
    muon_centrality_j1_j2 = FloatCol()
    met_phi_centrality = FloatCol()
    tau_j1_j2_phi_centrality = FloatCol()

    neff_pt = FloatCol()
    mass_all_jets = FloatCol()

    leadJetPt = FloatCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    weight = FloatCol()

    nvtx = IntCol()


    
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
    

class RecoTauMuBlock((RecoTau).prefix('tau_') + (RecoMuon).prefix('muon_')):

    @classmethod
    def set(cls, event, tree, tau, muon):
        """
        Misc variables
        """
        tree.mass_vis_tau_muon = utils.Mvis(tau.Et, tau.phi, muon.pt, muon.phi)
        tree.mass2_vis_tau_muon = (tau.fourvect + muon.fourvect).M()
        tree.theta_tau_muon = tau.fourvect.Vect().Angle(muon.fourvect.Vect())
        tree.cos_theta_tau_muon = math.cos(tree.theta_tau_muon)
        tree.charge_product_tau_muon = tau.charge * muon.charge
        tree.pt_ratio_tau_muon = muon.fourvect.Pt()/tau.fourvect.Pt()

        #Set tau variables
        setattr(tree, 'tau_BDTJetScore', tau.BDTJetScore)
        setattr(tree, 'tau_JetBDTSigLoose', tau.JetBDTSigLoose)
        setattr(tree, 'tau_JetBDTSigMedium', tau.JetBDTSigMedium)
        setattr(tree, 'tau_JetBDTSigTight', tau.JetBDTSigTight)
        setattr(tree, 'tau_numTrack', tau.numTrack)
        setattr(tree, 'tau_charge', tau.charge)
        setattr(tree, 'tau_chpiemeovercaloeme', tau.calcVars_ChPiEMEOverCaloEME)
        getattr(tree, 'tau_fourvect').set_from(tau.fourvect)

        #Set muon variables
        setattr(tree, 'muon_charge', muon.charge)
        getattr(tree, 'muon_fourvect').set_from(muon.fourvect)

        #Calculate muon isolation
        muon_pt = muon.fourvect.Pt()
        track_iso = muon.ptcone40/muon_pt <= 0.06
        calo_iso  = muon.etcone20/muon_pt <= 0.04
        setattr(tree, 'muon_isolated', (track_iso and calo_iso))


class RecoTauElectronBlock((RecoTau).prefix('tau_') + (RecoElectron).prefix('el_')):

    @classmethod
    def set(cls, event, tree, tau, electron):
        """
        Misc variables
        """
        tree.mass_vis_tau_electron = utils.Mvis(tau.Et, tau.phi, electron.fourvect.Pt(), electron.fourvect.Phi())
        tree.mass2_vis_tau_electron = (tau.fourvect + electron.fourvect).M()
        tree.theta_tau_electron = tau.fourvect.Vect().Angle(electron.fourvect.Vect())
        tree.cos_theta_tau_electron = math.cos(tree.theta_tau_electron)
        tree.charge_product_tau_electron = tau.charge * electron.charge
        if tau.fourvect.Pt() > 0:
            tree.pt_ratio_tau_electron = electron.fourvect.Pt() / tau.fourvect.Pt()
        else:
            tree.pt_ratio_tau_electron = -1111.

        # Set tau variables
        setattr(tree, 'tau_BDTJetScore', tau.BDTJetScore)
        setattr(tree, 'tau_JetBDTSigLoose', tau.JetBDTSigLoose)
        setattr(tree, 'tau_JetBDTSigMedium', tau.JetBDTSigMedium)
        setattr(tree, 'tau_JetBDTSigTight', tau.JetBDTSigTight)
        setattr(tree, 'tau_numTrack', tau.numTrack)
        setattr(tree, 'tau_charge', tau.charge)
        setattr(tree, 'tau_chpiemeovercaloeme', tau.calcVars_ChPiEMEOverCaloEME)
        getattr(tree, 'tau_fourvect').set_from(tau.fourvect)

        # Set electron variables
        setattr(tree, 'el_charge', electron.charge)
        getattr(tree, 'el_fourvect').set_from(electron.fourvect)

        # Calculate electron isolation
        track_iso = ( electron.ptcone40 / electron.cl_et < 0.06 )
        calo_iso  = ( electron.Etcone20 / electron.cl_et < 0.08 )
        setattr(tree, 'el_isolated', (track_iso and calo_iso))



class RecoTauLepBlock((RecoTau).prefix('tau_') + (RecoLepton).prefix('lep_')):

    @classmethod
    def set(cls, event, tree, tau, lep, leptype, isMC = True):
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
            calo_iso  = lep.etcone20/muon_pt <= 0.04
            setattr(tree, 'lep_isolated', (track_iso and calo_iso))

        if leptype == 1: # Is an electron
            track_iso = ( lep.ptcone40 / lep.cl_et < 0.06 )
            calo_iso  = ( lep.Etcone20 / lep.cl_et < 0.08 )
            setattr(tree, 'lep_isolated', (track_iso and calo_iso))

        #Find out if the tau is matched to an electron
        TauIsEl = False
        
        if tau.numTrack == 1 and tau.fourvect.Pt() > 20*GeV and isMC:
            nMC = event.mc_n
            for i in range(0, nMC):
                if abs(event.mc_pdgId[i]) == 11 and event.mc_pt[i] > 8*GeV:
                    if utils.dR(event.mc_eta[i], event.mc_phi[i], tau.eta, tau.phi) < 0.2: TauIsEl = True

        setattr(tree, "tau_isTruthEl", TauIsEl)

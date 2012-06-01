"""
This module defines the output branches in the final ntuple
as TreeModels.
"""


from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy.types import *

from atlastools.utils import et2pt
from atlastools import utils

from ..models import MatchedObject, TrueTau

import math
import ROOT


class RecoTau(TreeModel):

    BDTJetScore = FloatCol()
    chpiemeovercaloeme = FloatCol()

    JetBDTSigLoose = BoolCol()
    JetBDTSigMedium = BoolCol()
    JetBDTSigTight = BoolCol()

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
    mass_vis_tau_electron = FloatCol()
    mass2_vis_tau_electron = FloatCol()
    theta_tau_electron = FloatCol()
    cos_theta_tau_electron = FloatCol()
    pt_ratio_tau_electron = FloatCol()
    dphi_met_electron = FloatCol()

    numJets = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    numVertices = IntCol()
    HT = FloatCol()
    cutflow = IntCol()

    charge_product_tau_electron = IntCol()
    mass_transverse_met_electron = FloatCol()
    mass_transverse_met_tau = FloatCol()
    ddr_tau_electron = FloatCol()

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
    

class EventVariables(TreeModel):

    mass_collinear_tau_muon = FloatCol()
    tau_x = FloatCol()
    muon_x = FloatCol()
    mass_mmc_tau_muon = FloatCol()
    mass_vis_tau_muon = FloatCol()
    mass2_vis_tau_muon = FloatCol()
    theta_tau_muon = FloatCol()
    cos_theta_tau_muon = FloatCol()
    pt_ratio_tau_muon = FloatCol()
    dphi_met_muon = FloatCol()

    numJets = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    numVertices = IntCol()
    HT = FloatCol()
    cutflow = IntCol()

    charge_product_tau_muon = IntCol()
    mass_transverse_met_muon = FloatCol()
    mass_transverse_met_tau = FloatCol()
    ddr_tau_muon = FloatCol()

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

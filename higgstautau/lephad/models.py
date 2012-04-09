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


class RecoMET(TreeModel):

    MET_vect = Vector2
    MET_sig  = FloatCol()
    MET      = FloatCol()


class RecoJet(TreeModel):

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector

    jvtxf = FloatCol()


class EventVariables(TreeModel):

    mass_collinear_tau_muon = FloatCol()
    tau_x = FloatCol()
    muon_x = FloatCol()
    mass_mmc_tau_muon = FloatCol()
    mass_vis_tau_muon = FloatCol()
    mass2_vis_tau_muon = FloatCol()
    theta_tau_muon = FloatCol()
    cos_theta_tau_muon = FloatCol()

    numJets = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    numVertices = IntCol()
    HT = FloatCol()
    cutflow = IntCol()

    charge_product_tau_muon = IntCol()
    mass_transverse_met_muon = FloatCol()
    mass_j1_j2 = FloatCol()
    eta_product_j1_j2 = FloatCol()
    eta_delta_j1_j2 = FloatCol()
    tau_centrality_j1_j2 = FloatCol()
    muon_centrality_j1_j2 = FloatCol()
    met_phi_centrality = FloatCol()

    neff_pt = FloatCol()
    mass_all_jets = FloatCol()

    leadJetPt = FloatCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    weight = FloatCol()
    


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

        #Set tau variables
        setattr(tree, 'tau_BDTJetScore', tau.BDTJetScore)
        setattr(tree, 'tau_JetBDTSigLoose', tau.JetBDTSigLoose)
        setattr(tree, 'tau_JetBDTSigMedium', tau.JetBDTSigMedium)
        setattr(tree, 'tau_JetBDTSigTight', tau.JetBDTSigTight)
        setattr(tree, 'tau_numTrack', tau.numTrack)
        setattr(tree, 'tau_charge', tau.charge)
        getattr(tree, 'tau_fourvect').set_from(tau.fourvect)

        #Set muon variables
        setattr(tree, 'muon_charge', muon.charge)
        getattr(tree, 'muon_fourvect').set_from(muon.fourvect)

        #Calculate muon isolation
        muon_pt = muon.fourvect.Pt()
        track_iso = muon.ptcone40/muon_pt <= 0.06
        calo_iso  = muon.etcone20/muon_pt <= 0.04
        setattr(tree, 'muon_isolated', (track_iso and calo_iso))

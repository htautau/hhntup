"""
This module defines the output branches in the final ntuple
as TreeModels.
"""


from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzRotation, \
        LorentzVector, Vector3, Vector2
from rootpy.types import *

from atlastools.utils import et2pt
from atlastools import utils

from ..models import MatchedObject, TrueTau, FourMomentum

import math
import ROOT


class EventVariables(TreeModel):

    # event weight given by the PileupReweighting tool
    pileup_weight = FloatCol(default=1.)
    mc_weight = FloatCol(default=1.)
    ggf_weight = FloatCol(default=1.)

    # true if both taus pass ID requirements
    taus_pass = BoolCol()

    tau_trigger_match_error = BoolCol(default=False)

    theta_tau1_tau2 = FloatCol()
    cos_theta_tau1_tau2 = FloatCol()
    tau1_x = FloatCol()
    tau2_x = FloatCol()
    tau_x_product = FloatCol()
    tau_x_sum = FloatCol()
    tau_pt_ratio = FloatCol()
    tau_centrality_product = FloatCol()

    MET_centrality = FloatCol()
    MET_centrality_boosted = FloatCol()

    mass_collinear_tau1_tau2 = FloatCol()

    # new ditaumass package mass
    mass_dtm_tau1_tau2 = DoubleCol()
    mass_dtm_tau1_tau2_scan = DoubleCol()

    mass2_vis_tau1_tau2 = FloatCol()
    mass_vis_tau1_tau2 = FloatCol()

    mass_jet1_jet2 = FloatCol()
    mass_vis_true_tau1_tau2 = FloatCol()
    mass_true_quark1_quark2 = FloatCol()

    # did both taus come from the same vertex?
    tau_same_vertex = BoolCol()

    dR_quarks = FloatCol()
    dR_truetaus = FloatCol()
    dR_taus = FloatCol()
    dR_jets = FloatCol()
    dR_quark_tau = FloatCol()
    dR_tau1_tau2 = FloatCol()
    dEta_tau1_tau2 = FloatCol()
    dPhi_tau1_tau2 = FloatCol()

    dEta_quarks = FloatCol()
    dEta_jets = FloatCol()
    dEta_jets_boosted = FloatCol()
    eta_product_jets = FloatCol()
    eta_product_jets_boosted = FloatCol()

    numJets25 = IntCol()
    numJets = IntCol()
    nonisolatedjet = BoolCol()

    MET = FloatCol()
    MET_x = FloatCol()
    MET_y = FloatCol()
    MET_phi = FloatCol()
    MET_sig = FloatCol()
    MET_vec = Vector2
    dPhi_tau1_MET = FloatCol()
    dPhi_tau2_MET = FloatCol()
    dPhi_min_tau_MET = FloatCol()
    MET_bisecting = BoolCol()
    sumET = FloatCol()

    # MMC mass for all methods
    mmc0_mass = FloatCol()
    mmc1_mass = FloatCol()
    mmc2_mass = FloatCol()

    mmc0_MET = FloatCol()
    mmc0_MET_x = FloatCol()
    mmc0_MET_y = FloatCol()
    mmc0_MET_phi = FloatCol()
    mmc0_MET_vec = Vector2
    mmc1_MET = FloatCol()
    mmc1_MET_x = FloatCol()
    mmc1_MET_y = FloatCol()
    mmc1_MET_phi = FloatCol()
    mmc1_MET_vec = Vector2
    mmc2_MET = FloatCol()
    mmc2_MET_x = FloatCol()
    mmc2_MET_y = FloatCol()
    mmc2_MET_phi = FloatCol()
    mmc2_MET_vec = Vector2

    mmc0_resonance = LorentzVector
    mmc0_resonance_pt = FloatCol()
    mmc1_resonance = LorentzVector
    mmc1_resonance_pt = FloatCol()
    mmc2_resonance = LorentzVector
    mmc2_resonance_pt = FloatCol()

    jet_transformation = LorentzRotation
    jet_beta = Vector3
    parton_beta = Vector3

    error = BoolCol()
    cutflow = IntCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    sphericity_boosted = FloatCol()
    aplanarity_boosted = FloatCol()

    sum_pt = FloatCol()
    sum_pt_full = FloatCol()
    vector_sum_pt = FloatCol()

    ntrack_pv = IntCol()
    ntrack_nontau_pv = IntCol()


class EmbeddingBlock(TreeModel):

    embedding_isolation = IntCol()


class RecoTau(FourMomentum):

    BDTJetScore = FloatCol()
    BDTEleScore = FloatCol()

    JetBDTSigLoose = BoolCol()
    JetBDTSigMedium = BoolCol()
    JetBDTSigTight = BoolCol()

    nPi0 = IntCol()
    seedCalo_numTrack = IntCol()
    numTrack = IntCol()
    charge = IntCol()
    jvtxf = FloatCol()
    seedCalo_centFrac = FloatCol()

    centrality = FloatCol()
    centrality_boosted = FloatCol()

    # efficiency scale factor if matches truth
    efficiency_scale_factor = FloatCol(default=1.)
    efficiency_scale_factor_high = FloatCol(default=1.)
    efficiency_scale_factor_low = FloatCol(default=1.)

    # fake rate scale factor for taus that do not match truth
    fakerate_scale_factor = FloatCol(default=1.)
    fakerate_scale_factor_high = FloatCol(default=1.)
    fakerate_scale_factor_low = FloatCol(default=1.)

    # fake rate reco scale factor for taus that do not match truth
    fakerate_scale_factor_reco = FloatCol(default=1.)
    fakerate_scale_factor_reco_high = FloatCol(default=1.)
    fakerate_scale_factor_reco_low = FloatCol(default=1.)

    # trigger efficiency correction
    trigger_scale_factor = FloatCol(default=1.)
    trigger_scale_factor_high = FloatCol(default=1.)
    trigger_scale_factor_low = FloatCol(default=1.)

    trigger_match_thresh = IntCol(default=0)

    # overlap checking
    min_dr_jet = FloatCol(default=9999)

    # track recounting
    numTrack_recounted = IntCol(default=-1)

    # vertex association
    vertex_prob = FloatCol()


class RecoJet(FourMomentum):

    jvtxf = FloatCol()
    BDTJetScore = FloatCol()


class RecoTauBlock((RecoTau + MatchedObject).prefix('tau1_') +
                   (RecoTau + MatchedObject).prefix('tau2_')):

    @classmethod
    def set(cls, event, tree, tau1, tau2):

        tree.theta_tau1_tau2 = abs(tau1.fourvect.Angle(tau2.fourvect))
        tree.cos_theta_tau1_tau2 = math.cos(tree.theta_tau1_tau2)
        tree.dR_tau1_tau2 = tau1.fourvect.DeltaR(tau2.fourvect)
        tree.dEta_tau1_tau2 = abs(tau2.eta - tau1.eta)
        # leading pt over subleading pt
        tree.tau_pt_ratio = tau1.pt / tau2.pt
        tree.tau_centrality_product = tau1.centrality * tau2.centrality

        for outtau, intau in [(tree.tau1, tau1), (tree.tau2, tau2)]:

            FourMomentum.set(outtau, intau)

            outtau.BDTJetScore = intau.BDTJetScore
            outtau.BDTEleScore = intau.BDTEleScore

            outtau.JetBDTSigLoose = intau.JetBDTSigLoose
            outtau.JetBDTSigMedium = intau.JetBDTSigMedium
            outtau.JetBDTSigTight = intau.JetBDTSigTight

            outtau.nPi0 = intau.nPi0
            outtau.seedCalo_numTrack = intau.seedCalo_numTrack
            outtau.numTrack = intau.numTrack
            outtau.charge = intau.charge
            outtau.jvtxf = intau.jet_jvtxf
            outtau.seedCalo_centFrac = intau.seedCalo_centFrac

            outtau.centrality = intau.centrality
            outtau.centrality = intau.centrality_boosted

            if intau.matched:
                outtau.efficiency_scale_factor = intau.efficiency_scale_factor
                outtau.efficiency_scale_factor_high = intau.efficiency_scale_factor_high
                outtau.efficiency_scale_factor_low = intau.efficiency_scale_factor_low

                outtau.trigger_scale_factor = intau.trigger_scale_factor
                outtau.trigger_scale_factor_high = intau.trigger_scale_factor_high
                outtau.trigger_scale_factor_low = intau.trigger_scale_factor_low

            else:
                outtau.fakerate_scale_factor = intau.fakerate_scale_factor
                outtau.fakerate_scale_factor_high = intau.fakerate_scale_factor_high
                outtau.fakerate_scale_factor_low = intau.fakerate_scale_factor_low

                outtau.fakerate_scale_factor_reco = intau.fakerate_scale_factor_reco
                outtau.fakerate_scale_factor_reco_high = intau.fakerate_scale_factor_reco_high
                outtau.fakerate_scale_factor_reco_low = intau.fakerate_scale_factor_reco_low

            outtau.trigger_match_thresh = intau.trigger_match_thresh

            outtau.matched = intau.matched
            outtau.matched_dR = intau.matched_dR
            outtau.matched_collision = intau.matched_collision
            outtau.min_dr_jet = intau.min_dr_jet

            # track recounting
            outtau.numTrack_recounted = intau.numTrack_recounted

            # tau vertex association
            outtau.vertex_prob = intau.vertex_prob


class RecoJetBlock((RecoJet + MatchedObject).prefix('jet1_') +
                   (RecoJet + MatchedObject).prefix('jet2_')):

    @classmethod
    def set(cls, tree, jet1, jet2=None):

        FourMomentum.set(tree.jet1, jet1)
        tree.jet1_jvtxf = jet1.jvtxf

        if jet2 is not None:
            FourMomentum.set(tree.jet2, jet2)
            tree.jet2_jvtxf = jet2.jvtxf

            tree.mass_jet1_jet2 = (jet1.fourvect + jet2.fourvect).M()

            tree.dEta_jets = abs(
                    jet1.fourvect.Eta() - jet2.fourvect.Eta())
            tree.dEta_jets_boosted = abs(
                    jet1.fourvect_boosted.Eta() - jet2.fourvect_boosted.Eta())

            tree.eta_product_jets = jet1.fourvect.Eta() * jet2.fourvect.Eta()
            tree.eta_product_jets_boosted = (jet1.fourvect_boosted.Eta() *
                                             jet2.fourvect_boosted.Eta())


class TrueTauBlock((TrueTau + MatchedObject).prefix('trueTau1_') +
                   (TrueTau + MatchedObject).prefix('trueTau2_')):

    @classmethod
    def set(cls, tree, index, tau):

        setattr(tree, 'trueTau%i_nProng' % index, tau.nProng)
        setattr(tree, 'trueTau%i_nPi0' % index, tau.nPi0)
        setattr(tree, 'trueTau%i_charge' % index, tau.charge)

        fourvect = getattr(tree, 'trueTau%i_fourvect' % index)
        fourvect.SetPtEtaPhiM(
            tau.pt,
            tau.eta,
            tau.phi,
            tau.m)

        fourvect_boosted = getattr(tree, 'trueTau%i_fourvect_boosted' % index)
        fourvect_boosted.set_from(fourvect)
        fourvect_boosted.Boost(tree.parton_beta * -1)

        fourvect_vis = getattr(tree, 'trueTau%i_fourvect_vis' % index)
        try:
            fourvect_vis.SetPtEtaPhiM(
                et2pt(tau.vis_Et, tau.vis_eta, tau.vis_m),
                tau.vis_eta,
                tau.vis_phi,
                tau.vis_m)
        except ValueError:
            print "DOMAIN ERROR ON TRUTH 4VECT"
            print tau.vis_Et, tau.vis_eta, tau.vis_m
        else:
            fourvect_vis_boosted = getattr(tree, 'trueTau%i_fourvect_vis_boosted' % index)
            fourvect_vis_boosted.set_from(fourvect_vis)
            fourvect_vis_boosted.Boost(tree.parton_beta * -1)

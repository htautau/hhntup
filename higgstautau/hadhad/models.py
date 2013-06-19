"""
This module defines the output branches in the final ntuple
as TreeModels.
"""

from rootpy.tree import TreeModel, FloatCol, IntCol, DoubleCol, BoolCol
from rootpy.math.physics.vector import (
    LorentzRotation, LorentzVector, Vector3, Vector2)
from rootpy import log
ignore_warning = log['/ROOT.TVector3.PseudoRapidity'].ignore(
        '.*transvers momentum.*')

from atlastools.utils import et2pt
from atlastools import datasets
from atlastools import utils

from ..models import MatchedObject

import math
import ROOT
from ROOT import TLorentzVector


class FourMomentum(TreeModel):

    pt = FloatCol()
    p = FloatCol()
    et = FloatCol()
    e = FloatCol()
    eta = FloatCol(default=-1111)
    phi = FloatCol(default=-1111)
    m = FloatCol()

    @classmethod
    def set(cls, this, other):

        if isinstance(other, TLorentzVector):
            vect = other
        else:
            vect = other.fourvect
        this.pt = vect.Pt()
        this.p = vect.P()
        this.et = vect.Et()
        this.e = vect.E()
        this.m = vect.M()
        with ignore_warning:
            this.phi = vect.Phi()
            this.eta = vect.Eta()


class TrueTau(TreeModel):

    nProng = IntCol(default=-1111)
    nPi0 = IntCol(default=-1111)
    charge = IntCol()


class MMCOutput(FourMomentum.prefix('resonance_')):

    mass = FloatCol()
    MET = FloatCol()
    MET_x = FloatCol()
    MET_y = FloatCol()
    MET_phi = FloatCol()
    MET_vec = Vector2


class MMCModel(MMCOutput.prefix('mmc0_'),
               MMCOutput.prefix('mmc1_'),
               MMCOutput.prefix('mmc2_')):
    pass


class MassModel(MMCModel):

    mass_collinear_tau1_tau2 = FloatCol()
    mass_vis_tau1_tau2 = FloatCol()
    mass_jet1_jet2 = FloatCol()
    mass_vis_true_tau1_tau2 = FloatCol()
    mass_true_quark1_quark2 = FloatCol()


class METModel(TreeModel):

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
    MET_centrality = FloatCol()
    #MET_centrality_boosted = FloatCol()


class EmbeddingModel(TreeModel):

    # embedding corrections
    embedding_reco_unfold = FloatCol(default=1.)
    embedding_dimuon_mass = FloatCol()
    embedding_isolation = IntCol()


class RecoTau(FourMomentum):

    index = IntCol(default=-1)

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
    #centrality_boosted = FloatCol()

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

    collinear_momentum_fraction = FloatCol()


class RecoJet(FourMomentum):

    index = IntCol(default=-1)

    jvtxf = FloatCol()
    BDTJetScore = FloatCol()


class RecoTauBlock((RecoTau + MatchedObject).prefix('tau1_') +
                   (RecoTau + MatchedObject).prefix('tau2_')):

    # true if both taus pass ID requirements
    taus_pass = BoolCol()

    tau_trigger_match_error = BoolCol(default=False)

    # did both taus come from the same vertex?
    tau_same_vertex = BoolCol()

    theta_tau1_tau2 = FloatCol()
    cos_theta_tau1_tau2 = FloatCol()
    tau1_x = FloatCol()
    tau2_x = FloatCol()

    tau_x_product = FloatCol()
    tau_x_sum = FloatCol()
    tau_pt_ratio = FloatCol()
    tau_centrality_product = FloatCol()

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

            outtau.index = intau.index

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
            #outtau.centrality_boosted = intau.centrality_boosted

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

            outtau.collinear_momentum_fraction = intau.collinear_momentum_fraction


class RecoJetBlock((RecoJet + MatchedObject).prefix('jet1_') +
                   (RecoJet + MatchedObject).prefix('jet2_')):

    #jet_transformation = LorentzRotation
    #jet_beta = Vector3
    #parton_beta = Vector3
    numJets = IntCol()
    nonisolatedjet = BoolCol()

    @classmethod
    def set(cls, tree, jet1, jet2=None):

        FourMomentum.set(tree.jet1, jet1)
        tree.jet1_jvtxf = jet1.jvtxf
        tree.jet1_index = jet1.index

        if jet2 is not None:

            FourMomentum.set(tree.jet2, jet2)
            tree.jet2_jvtxf = jet2.jvtxf
            tree.jet2_index = jet2.index

            tree.mass_jet1_jet2 = (jet1.fourvect + jet2.fourvect).M()

            tree.dEta_jets = abs(
                    jet1.fourvect.Eta() - jet2.fourvect.Eta())
            #tree.dEta_jets_boosted = abs(
            #        jet1.fourvect_boosted.Eta() - jet2.fourvect_boosted.Eta())

            tree.eta_product_jets = jet1.fourvect.Eta() * jet2.fourvect.Eta()
            #tree.eta_product_jets_boosted = (jet1.fourvect_boosted.Eta() *
            #                                 jet2.fourvect_boosted.Eta())


class TrueTauBlock((TrueTau + MatchedObject).prefix('truetau1_') +
                   (TrueTau + MatchedObject).prefix('truetau2_')):

    @classmethod
    def set(cls, tree, index, tau):

        setattr(tree, 'truetau%i_nProng' % index, tau.nProng)
        setattr(tree, 'truetau%i_nPi0' % index, tau.nPi0)
        setattr(tree, 'truetau%i_charge' % index, tau.charge)

        fourvect = getattr(tree, 'truetau%i_fourvect' % index)
        fourvect.SetPtEtaPhiM(
            tau.pt,
            tau.eta,
            tau.phi,
            tau.m)

        fourvect_boosted = getattr(tree, 'truetau%i_fourvect_boosted' % index)
        fourvect_boosted.set_from(fourvect)
        fourvect_boosted.Boost(tree.parton_beta * -1)

        fourvect_vis = getattr(tree, 'truetau%i_fourvect_vis' % index)
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
            fourvect_vis_boosted = getattr(tree, 'truetau%i_fourvect_vis_boosted' % index)
            fourvect_vis_boosted.set_from(fourvect_vis)
            fourvect_vis_boosted.Boost(tree.parton_beta * -1)


class EventModel(TreeModel):

    RunNumber = IntCol()
    number_of_good_vertices = IntCol()
    averageIntPerXing_scaled = FloatCol()
    actualIntPerXing_scaled = FloatCol()

    # event weight given by the PileupReweighting tool
    pileup_weight = FloatCol(default=1.)
    mc_weight = FloatCol(default=1.)
    ggf_weight = FloatCol(default=1.)

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
    #dEta_jets_boosted = FloatCol()
    eta_product_jets = FloatCol()
    #eta_product_jets_boosted = FloatCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    #sphericity_boosted = FloatCol()
    #aplanarity_boosted = FloatCol()

    sum_pt = FloatCol()
    sum_pt_full = FloatCol()
    vector_sum_pt = FloatCol()
    resonance_pt = FloatCol()

    ntrack_pv = IntCol()
    ntrack_nontau_pv = IntCol()

    error = BoolCol()


def get_model(datatype, name, prefix=None):

    model = EventModel + MassModel + METModel + RecoTauBlock + RecoJetBlock
    #if datatype in (datasets.MC, datasets.EMBED):
    #    model += TrueTauBlock
    if datatype == datasets.EMBED:
        model += EmbeddingModel
    #if datatype == datasets.MC and 'VBF' in name:
    #    # add branches for VBF Higgs associated partons
    #    model += PartonBlock
    if prefix is not None:
        return model.prefix(prefix)
    return model

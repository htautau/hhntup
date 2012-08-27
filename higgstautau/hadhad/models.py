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
    BDTEleScore = FloatCol()

    JetBDTSigLoose = BoolCol()
    JetBDTSigMedium = BoolCol()
    JetBDTSigTight = BoolCol()

    nPi0 = IntCol()
    seedCalo_numTrack = IntCol()
    numTrack = IntCol()
    charge = IntCol()

    centrality = FloatCol()
    centrality_boosted = FloatCol()

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector

    # efficiency scale factor if matches truth
    efficiency_scale_factor = FloatCol(default=1.)
    efficiency_scale_factor_high = FloatCol(default=1.)
    efficiency_scale_factor_low = FloatCol(default=1.)

    # fake rate scale factor for taus that do not match truth
    fakerate_scale_factor = FloatCol(default=1.)
    fakerate_scale_factor_high = FloatCol(default=1.)
    fakerate_scale_factor_low = FloatCol(default=1.)

    # trigger efficiency correction
    trigger_scale_factor = FloatCol(default=1.)
    trigger_scale_factor_high = FloatCol(default=1.)
    trigger_scale_factor_low = FloatCol(default=1.)


class EventVariables(TreeModel):

    # the category (2jet or 01jet)
    category = IntCol()

    # event weight given by the PileupReweighting tool
    pileup_weight = FloatCol(default=1.)
    mc_weight = FloatCol(default=1.)
    ggf_weight = FloatCol(default=1.)

    theta_tau1_tau2 = FloatCol()
    cos_theta_tau1_tau2 = FloatCol()
    tau1_x = FloatCol()
    tau2_x = FloatCol()

    MET_centrality = FloatCol()
    MET_centrality_boosted = FloatCol()

    mass_collinear_tau1_tau2 = FloatCol()
    # MMC mass
    mass_mmc_tau1_tau2 = FloatCol()

    # new ditaumass package mass
    mass_dtm_tau1_tau2 = DoubleCol()
    mass_dtm_tau1_tau2_scan = DoubleCol()

    mass2_vis_tau1_tau2 = FloatCol()
    mass_vis_tau1_tau2 = FloatCol()

    mass_jet1_jet2 = FloatCol()
    mass_vis_true_tau1_tau2 = FloatCol()
    mass_true_quark1_quark2 = FloatCol()

    dR_quarks = FloatCol()
    dR_truetaus = FloatCol()
    dR_taus = FloatCol()
    dR_jets = FloatCol()
    dR_quark_tau = FloatCol()
    dR_tau1_tau2 = FloatCol()
    dPhi_tau1_tau2 = FloatCol()

    dEta_quarks = FloatCol()
    dEta_jets = FloatCol()
    dEta_jets_boosted = FloatCol()
    eta_product_jets = FloatCol()
    eta_product_jets_boosted = FloatCol()

    numJets = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    MET = FloatCol()
    MET_phi = FloatCol()
    MET_mmc = FloatCol()
    HT = FloatCol()
    MET_sig = FloatCol()
    error = BoolCol()
    jet_transformation = LorentzRotation
    jet_beta = Vector3
    parton_beta = Vector3
    cutflow = IntCol()

    sphericity = FloatCol()
    aplanarity = FloatCol()

    sphericity_full = FloatCol()
    aplanarity_full = FloatCol()

    sphericity_boosted = FloatCol()
    aplanarity_boosted = FloatCol()

    higgs_pt = FloatCol()
    sum_pt = FloatCol()
    sum_pt_full = FloatCol()


class RecoJet(TreeModel):

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector

    jvtxf = FloatCol()
    BDTJetScore = FloatCol()


class RecoTauBlock((RecoTau + MatchedObject).prefix('tau1_') + (RecoTau + MatchedObject).prefix('tau2_')):

    @classmethod
    def set(cls, event, tree, tau1, tau2):

        if tau1 is not None and tau2 is not None:
            tree.mass_vis_tau1_tau2 = utils.Mvis(tau1.Et, tau1.seedCalo_phi, tau2.Et, tau2.seedCalo_phi)
            tree.mass2_vis_tau1_tau2 = (tau1.fourvect + tau2.fourvect).M()
            tree.theta_tau1_tau2 = tau1.fourvect.Vect().Angle(tau2.fourvect.Vect())
            tree.cos_theta_tau1_tau2 = math.cos(tree.theta_tau1_tau2)
            tree.dR_tau1_tau2 = tau1.fourvect.DeltaR(tau2.fourvect)
            tree.dPhi_tau1_tau2 = abs(tau1.fourvect.DeltaPhi(tau2.fourvect))

        for i, tau in zip((1,2), (tau1, tau2)):
            if tau is None:
                continue
            fourvect = tau.fourvect
            setattr(tree, 'tau%i_BDTJetScore' % i, tau.BDTJetScore)
            setattr(tree, 'tau%i_BDTEleScore' % i, tau.BDTEleScore)

            setattr(tree, 'tau%i_JetBDTSigLoose' % i, tau.JetBDTSigLoose)
            setattr(tree, 'tau%i_JetBDTSigMedium' % i, tau.JetBDTSigMedium)
            setattr(tree, 'tau%i_JetBDTSigTight' % i, tau.JetBDTSigTight)

            setattr(tree, 'tau%i_nPi0' % i, tau.nPi0)
            setattr(tree, 'tau%i_seedCalo_numTrack' % i, tau.seedCalo_numTrack)
            setattr(tree, 'tau%i_numTrack' % i, tau.numTrack)
            setattr(tree, 'tau%i_charge' % i, tau.charge)
            getattr(tree, 'tau%i_fourvect' % i).set_from(tau.fourvect)
            tau.fourvect_boosted.set_from(tau.fourvect)
            tau.fourvect_boosted.Boost(tree.jet_beta * -1)
            getattr(tree, 'tau%i_fourvect_boosted' % i).set_from(tau.fourvect_boosted)

            setattr(tree, 'tau%i_centrality' % i, tau.centrality)
            setattr(tree, 'tau%i_centrality_boosted' % i, tau.centrality_boosted)

            setattr(tree, 'tau%i_efficiency_scale_factor' % i,
                    tau.efficiency_scale_factor)
            setattr(tree, 'tau%i_efficiency_scale_factor_high' % i,
                    tau.efficiency_scale_factor_high)
            setattr(tree, 'tau%i_efficiency_scale_factor_low' % i,
                    tau.efficiency_scale_factor_low)

            setattr(tree, 'tau%i_fakerate_scale_factor' % i,
                    tau.fakerate_scale_factor)
            setattr(tree, 'tau%i_fakerate_scale_factor_high' % i,
                    tau.fakerate_scale_factor_high)
            setattr(tree, 'tau%i_fakerate_scale_factor_low' % i,
                    tau.fakerate_scale_factor_low)

            setattr(tree, 'tau%i_trigger_scale_factor' % i,
                    tau.trigger_scale_factor)
            setattr(tree, 'tau%i_trigger_scale_factor_high' % i,
                    tau.trigger_scale_factor_high)
            setattr(tree, 'tau%i_trigger_scale_factor_low' % i,
                    tau.trigger_scale_factor_low)

            setattr(tree, 'tau%i_matched' % i, tau.matched)
            setattr(tree, 'tau%i_matched_dR' % i, tau.matched_dR)
            setattr(tree, 'tau%i_matched_collision' % i, tau.matched_collision)


class RecoJetBlock((RecoJet + MatchedObject).prefix('jet1_') + (RecoJet + MatchedObject).prefix('jet2_')):

    @classmethod
    def set(cls, tree, jet1, jet2=None):

        tree.jet1_fourvect.set_from(jet1.fourvect)
        tree.jet1_jvtxf = jet1.jvtxf

        if jet2 is not None:
            # the jets should already be sorted by eta
            # sort by eta
            # jet1, jet2 = sorted([jet1, jet2], key=lambda jet: jet.fourvect.Eta())

            # determine jet CoM frame
            beta = (jet1.fourvect + jet2.fourvect).BoostVector()
            tree.jet_beta.set_from(beta)

            jet1.fourvect_boosted.set_from(jet1.fourvect)
            jet2.fourvect_boosted.set_from(jet2.fourvect)
            jet1.fourvect_boosted.Boost(beta * -1)
            jet2.fourvect_boosted.Boost(beta * -1)

            # sort by transformed eta
            #jet1, jet2 = sorted([jet1, jet2], key=lambda jet: jet.fourvect_boosted.Eta())

            tree.mass_jet1_jet2 = (jet1.fourvect + jet2.fourvect).M()

            tree.jet2_fourvect.set_from(jet2.fourvect)

            tree.jet1_fourvect_boosted.set_from(jet1.fourvect_boosted)
            tree.jet2_fourvect_boosted.set_from(jet2.fourvect_boosted)

            tree.jet2_jvtxf = jet2.jvtxf

            tree.dEta_jets = abs(jet1.fourvect.Eta() - jet2.fourvect.Eta())
            tree.dEta_jets_boosted = abs(jet1.fourvect_boosted.Eta() - jet2.fourvect_boosted.Eta())

            tree.eta_product_jets = jet1.fourvect.Eta() * jet2.fourvect.Eta()
            tree.eta_product_jets_boosted = (jet1.fourvect_boosted.Eta() *
                                             jet2.fourvect_boosted.Eta())


class TrueTauBlock((TrueTau + MatchedObject).prefix('trueTau1_') + (TrueTau + MatchedObject).prefix('trueTau2_')):

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
        fourvect_vis.SetPtEtaPhiM(
            et2pt(tau.vis_Et, tau.vis_eta, tau.vis_m),
            tau.vis_eta,
            tau.vis_phi,
            tau.vis_m)

        fourvect_vis_boosted = getattr(tree, 'trueTau%i_fourvect_vis_boosted' % index)
        fourvect_vis_boosted.set_from(fourvect_vis)
        fourvect_vis_boosted.Boost(tree.parton_beta * -1)

"""
This module defines the output branches in the final ntuple
as TreeModels.
"""

from rootpy.tree import TreeModel, FloatCol, IntCol, DoubleCol, BoolCol
from rootpy.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy import stl

from .. import datasets
from .. import eventshapes
from ..models import FourMomentum, MatchedObject, MMCModel, TrueTau

import math


class MassModel(MMCModel):
    mass_collinear_tau1_tau2 = FloatCol()
    mass_vis_tau1_tau2 = FloatCol()
    mass_jet1_jet2 = FloatCol(default=-1E10)
    mass_tau1_tau2_jet1 = FloatCol(default=-1E10)
    #mass_vis_true_tau1_tau2 = FloatCol()
    #mass_true_quark1_quark2 = FloatCol()


class METModel(TreeModel):
    MET_et = FloatCol()
    MET_etx = FloatCol()
    MET_ety = FloatCol()
    MET_sumet = FloatCol()
    MET_phi = FloatCol()
    MET_sig = FloatCol()

    MET_et_original = FloatCol()
    MET_etx_original = FloatCol()
    MET_ety_original = FloatCol()
    MET_sumet_original = FloatCol()
    MET_phi_original = FloatCol()

    dPhi_tau1_MET = FloatCol()
    dPhi_tau2_MET = FloatCol()
    dPhi_min_tau_MET = FloatCol()
    MET_bisecting = BoolCol()

    MET_centrality = FloatCol(default=-1E10)
    MET_centrality_boosted = FloatCol(default=-1E10)


class EmbeddingModel(TreeModel):
    # embedding corrections
    embedding_reco_unfold = FloatCol(default=1.)
    embedding_trigger_weight = FloatCol(default=1.)
    embedding_dimuon_mass = FloatCol()
    embedding_isolation = IntCol()
    embedding_spin_weight = FloatCol(default=1.)


class RecoTau(FourMomentum):
    index = IntCol(default=-1)

    BDTJetScore = FloatCol()
    BDTEleScore = FloatCol()

    JetBDTSigLoose = BoolCol()
    JetBDTSigMedium = BoolCol()
    JetBDTSigTight = BoolCol()
    id = IntCol()

    nPi0 = IntCol()
    seedCalo_numTrack = IntCol()
    numTrack = IntCol()
    charge = IntCol()
    jvtxf = FloatCol()
    seedCalo_centFrac = FloatCol()

    BCHMedium = BoolCol()
    BCHTight = BoolCol()

    centrality = FloatCol(default=-1E10)
    centrality_boosted = FloatCol(default=-1E10)

    # efficiency scale factor if matches truth
    id_sf = FloatCol(default=1.)
    id_sf_high = FloatCol(default=1.)
    id_sf_low = FloatCol(default=1.)
    id_sf_stat_high = FloatCol(default=1.)
    id_sf_stat_low = FloatCol(default=1.)
    id_sf_sys_high = FloatCol(default=1.)
    id_sf_sys_low = FloatCol(default=1.)

    # trigger efficiency
    trigger_sf = FloatCol(default=1.)
    trigger_sf_high = FloatCol(default=1.)
    trigger_sf_low = FloatCol(default=1.)
    trigger_sf_mc_stat_high = FloatCol(default=1.)
    trigger_sf_mc_stat_low = FloatCol(default=1.)
    trigger_sf_data_stat_high = FloatCol(default=1.)
    trigger_sf_data_stat_low = FloatCol(default=1.)
    trigger_sf_stat_high = FloatCol(default=1.)
    trigger_sf_stat_low = FloatCol(default=1.)
    trigger_sf_stat_scale_high = FloatCol(default=1.)
    trigger_sf_stat_scale_low = FloatCol(default=1.)
    trigger_sf_sys_high = FloatCol(default=1.)
    trigger_sf_sys_low = FloatCol(default=1.)

    trigger_eff = FloatCol(default=1.)
    trigger_eff_high = FloatCol(default=1.)
    trigger_eff_low = FloatCol(default=1.)
    trigger_eff_stat_high = FloatCol(default=1.)
    trigger_eff_stat_low = FloatCol(default=1.)
    trigger_eff_stat_scale_high = FloatCol(default=1.)
    trigger_eff_stat_scale_low = FloatCol(default=1.)
    trigger_eff_sys_high = FloatCol(default=1.)
    trigger_eff_sys_low = FloatCol(default=1.)

    trigger_sf_stat_scale_PeriodA_high = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodA_low = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodBD_Barrel_high = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodBD_Barrel_low = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodBD_EndCap_high = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodBD_EndCap_low = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodEM_Barrel_high = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodEM_Barrel_low = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodEM_EndCap_high = FloatCol(default=1.)
    trigger_sf_stat_scale_PeriodEM_EndCap_low = FloatCol(default=1.)

    trigger_eff_stat_scale_PeriodA_high = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodA_low = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodBD_Barrel_high = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodBD_Barrel_low = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodBD_EndCap_high = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodBD_EndCap_low = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodEM_Barrel_high = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodEM_Barrel_low = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodEM_EndCap_high = FloatCol(default=1.)
    trigger_eff_stat_scale_PeriodEM_EndCap_low = FloatCol(default=1.)

    # fake rate scale factor for taus that do not match truth
    fakerate_sf = FloatCol(default=1.)
    fakerate_sf_high = FloatCol(default=1.)
    fakerate_sf_low = FloatCol(default=1.)

    # fake rate reco scale factor for taus that do not match truth
    fakerate_sf_reco = FloatCol(default=1.)
    fakerate_sf_reco_high = FloatCol(default=1.)
    fakerate_sf_reco_low = FloatCol(default=1.)

    #trigger_match_thresh = IntCol(default=0)

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
    BCHMedium = BoolCol()
    BCHTight = BoolCol()


class RecoTauBlock((RecoTau + MatchedObject).prefix('tau1_') +
                   (RecoTau + MatchedObject).prefix('tau2_')):

    #tau_trigger_match_error = BoolCol(default=False)

    # did both taus come from the same vertex?
    tau_same_vertex = BoolCol()

    dR_tau1_tau2 = FloatCol()
    dEta_tau1_tau2 = FloatCol()
    theta_tau1_tau2 = FloatCol()
    cos_theta_tau1_tau2 = FloatCol()
    tau_pt_ratio = FloatCol()
    # set in hhskim.py
    dPhi_tau1_tau2 = FloatCol()

    @classmethod
    def set(cls, event, tree, datatype, tau1, tau2, local=False):
        tree.theta_tau1_tau2 = abs(tau1.fourvect.Angle(tau2.fourvect))
        tree.cos_theta_tau1_tau2 = math.cos(tree.theta_tau1_tau2)
        tree.dR_tau1_tau2 = tau1.fourvect.DeltaR(tau2.fourvect)
        tree.dEta_tau1_tau2 = abs(tau2.eta - tau1.eta)
        # leading pt over subleading pt
        tree.tau_pt_ratio = tau1.pt / tau2.pt

        for outtau, intau in [(tree.tau1, tau1), (tree.tau2, tau2)]:

            outtau.index = intau.index
            outtau.id = intau.id

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
            #outtau.jvtxf = intau.jet_jvtxf
            outtau.seedCalo_centFrac = intau.seedCalo_centFrac

            outtau.centrality = intau.centrality
            outtau.centrality_boosted = intau.centrality_boosted

            if intau.matched:
                outtau.id_sf = intau.id_sf
                outtau.id_sf_high = intau.id_sf_high
                outtau.id_sf_low = intau.id_sf_low
                outtau.id_sf_stat_high = intau.id_sf_stat_high
                outtau.id_sf_stat_low = intau.id_sf_stat_low
                outtau.id_sf_sys_high = intau.id_sf_sys_high
                outtau.id_sf_sys_low = intau.id_sf_sys_low

                outtau.trigger_sf = intau.trigger_sf
                outtau.trigger_sf_high = intau.trigger_sf_high
                outtau.trigger_sf_low = intau.trigger_sf_low
                outtau.trigger_sf_mc_stat_high = intau.trigger_sf_mc_stat_high
                outtau.trigger_sf_mc_stat_low = intau.trigger_sf_mc_stat_low
                outtau.trigger_sf_data_stat_high = intau.trigger_sf_data_stat_high
                outtau.trigger_sf_data_stat_low = intau.trigger_sf_data_stat_low
                outtau.trigger_sf_stat_high = intau.trigger_sf_stat_high
                outtau.trigger_sf_stat_low = intau.trigger_sf_stat_low
                outtau.trigger_sf_stat_scale_high = intau.trigger_sf_stat_scale_high
                outtau.trigger_sf_stat_scale_low = intau.trigger_sf_stat_scale_low
                outtau.trigger_sf_sys_high = intau.trigger_sf_sys_high
                outtau.trigger_sf_sys_low = intau.trigger_sf_sys_low

                outtau.trigger_eff = intau.trigger_eff
                outtau.trigger_eff_high = intau.trigger_eff_high
                outtau.trigger_eff_low = intau.trigger_eff_low
                outtau.trigger_eff_stat_high = intau.trigger_eff_stat_high
                outtau.trigger_eff_stat_low = intau.trigger_eff_stat_low
                outtau.trigger_eff_stat_scale_high = intau.trigger_eff_stat_scale_high
                outtau.trigger_eff_stat_scale_low = intau.trigger_eff_stat_scale_low
                outtau.trigger_eff_sys_high = intau.trigger_eff_sys_high
                outtau.trigger_eff_sys_low = intau.trigger_eff_sys_low

                # partitioned trigger stat NP
                # same bins as in TrigTauEfficiency
                if 200804 <= tree.RunNumber <= 201556:
                    # period A eta inclusive
                    outtau.trigger_sf_stat_scale_PeriodA_high = intau.trigger_sf_stat_scale_high
                    outtau.trigger_sf_stat_scale_PeriodA_low = intau.trigger_sf_stat_scale_low
                    outtau.trigger_eff_stat_scale_PeriodA_high = intau.trigger_eff_stat_scale_high
                    outtau.trigger_eff_stat_scale_PeriodA_low = intau.trigger_eff_stat_scale_low
                elif 202660 <= tree.RunNumber <= 209025:
                    # period B-D
                    if abs(intau.eta) <= 1.5:
                        outtau.trigger_sf_stat_scale_PeriodBD_Barrel_high = intau.trigger_sf_stat_scale_high
                        outtau.trigger_sf_stat_scale_PeriodBD_Barrel_low = intau.trigger_sf_stat_scale_low
                        outtau.trigger_eff_stat_scale_PeriodBD_Barrel_high = intau.trigger_eff_stat_scale_high
                        outtau.trigger_eff_stat_scale_PeriodBD_Barrel_low = intau.trigger_eff_stat_scale_low
                    else:
                        outtau.trigger_sf_stat_scale_PeriodBD_EndCap_high = intau.trigger_sf_stat_scale_high
                        outtau.trigger_sf_stat_scale_PeriodBD_EndCap_low = intau.trigger_sf_stat_scale_low
                        outtau.trigger_eff_stat_scale_PeriodBD_EndCap_high = intau.trigger_eff_stat_scale_high
                        outtau.trigger_eff_stat_scale_PeriodBD_EndCap_low = intau.trigger_eff_stat_scale_low
                elif 209074 <= tree.RunNumber <= 216432:
                    # period E-M
                    if abs(intau.eta) <= 1.5:
                        outtau.trigger_sf_stat_scale_PeriodEM_Barrel_high = intau.trigger_sf_stat_scale_high
                        outtau.trigger_sf_stat_scale_PeriodEM_Barrel_low = intau.trigger_sf_stat_scale_low
                        outtau.trigger_eff_stat_scale_PeriodEM_Barrel_high = intau.trigger_eff_stat_scale_high
                        outtau.trigger_eff_stat_scale_PeriodEM_Barrel_low = intau.trigger_eff_stat_scale_low
                    else:
                        outtau.trigger_sf_stat_scale_PeriodEM_EndCap_high = intau.trigger_sf_stat_scale_high
                        outtau.trigger_sf_stat_scale_PeriodEM_EndCap_low = intau.trigger_sf_stat_scale_low
                        outtau.trigger_eff_stat_scale_PeriodEM_EndCap_high = intau.trigger_eff_stat_scale_high
                        outtau.trigger_eff_stat_scale_PeriodEM_EndCap_low = intau.trigger_eff_stat_scale_low

            else:
                outtau.fakerate_sf = intau.fakerate_sf
                outtau.fakerate_sf_high = intau.fakerate_sf_high
                outtau.fakerate_sf_low = intau.fakerate_sf_low

                outtau.fakerate_sf_reco = intau.fakerate_sf_reco
                outtau.fakerate_sf_reco_high = intau.fakerate_sf_reco_high
                outtau.fakerate_sf_reco_low = intau.fakerate_sf_reco_low

            #outtau.trigger_match_thresh = intau.trigger_match_thresh

            outtau.matched = intau.matched
            outtau.matched_dR = intau.matched_dR
            # outtau.matched_collision = intau.matched_collision
            outtau.min_dr_jet = intau.min_dr_jet

            if not local:
                # track recounting
                # the track branches are removed by the skim, so this should
                # only be set in the skim and cannot be recomputed on the skim
                outtau.numTrack_recounted = intau.numTrack_recounted
                # BCH cleaning only computed in skim
                outtau.BCHMedium = intau.BCHMedium
                outtau.BCHTight = intau.BCHTight

            # tau vertex association
            outtau.vertex_prob = intau.vertex_prob

            outtau.collinear_momentum_fraction = intau.collinear_momentum_fraction


class RecoJetBlock(RecoJet.prefix('jet1_') +
                   RecoJet.prefix('jet2_') +
                   RecoJet.prefix('jet3_')):

    dEta_jets = FloatCol(default=-1)
    dEta_jets_boosted = FloatCol(default=-1)
    eta_product_jets = FloatCol(default=-1E10)
    eta_product_jets_boosted = FloatCol(default=-1E10)
    #jet_transformation = LorentzRotation
    jet_beta = Vector3
    numJets = IntCol()
    nonisolatedjet = BoolCol()
    jet3_centrality = FloatCol(default=-1E10)
    jet3_centrality_boosted = FloatCol(default=-1E10)

    @classmethod
    def set(cls, tree, jet1, jet2=None, jet3=None, local=False):

        if jet1 is not None:
            FourMomentum.set(tree.jet1, jet1)
            tree.jet1_jvtxf = jet1.jvtxf
            tree.jet1_index = jet1.index

            if not local:
                # only computed in skim!
                tree.jet1_BCHMedium = jet1.BCHMedium
                tree.jet1_BCHTight = jet1.BCHTight

        elif local:
            # zero the fourvect
            # avoid ghost values from skim if jet selection changed
            tree['jet1_pt'].reset()
            tree['jet1_p'].reset()
            tree['jet1_et'].reset()
            tree['jet1_e'].reset()
            tree['jet1_m'].reset()
            tree['jet1_phi'].reset()
            tree['jet1_eta'].reset()

        if jet2 is not None:
            FourMomentum.set(tree.jet2, jet2)
            tree.jet2_jvtxf = jet2.jvtxf
            tree.jet2_index = jet2.index

            if not local:
                # only computed in skim!
                tree.jet2_BCHMedium = jet2.BCHMedium
                tree.jet2_BCHTight = jet2.BCHTight

            tree.mass_jet1_jet2 = (jet1.fourvect + jet2.fourvect).M()

            tree.dEta_jets = abs(
                jet1.fourvect.Eta() - jet2.fourvect.Eta())
            tree.dEta_jets_boosted = abs(
                jet1.fourvect_boosted.Eta() - jet2.fourvect_boosted.Eta())

            tree.eta_product_jets = jet1.fourvect.Eta() * jet2.fourvect.Eta()
            tree.eta_product_jets_boosted = (
                jet1.fourvect_boosted.Eta() * jet2.fourvect_boosted.Eta())
        elif local:
            # zero the fourvect
            # avoid ghost values from skim if jet selection changed
            tree['jet2_pt'].reset()
            tree['jet2_p'].reset()
            tree['jet2_et'].reset()
            tree['jet2_e'].reset()
            tree['jet2_m'].reset()
            tree['jet2_phi'].reset()
            tree['jet2_eta'].reset()

        if jet3 is not None:
            FourMomentum.set(tree.jet3, jet3)
            tree.jet3_jvtxf = jet3.jvtxf
            tree.jet3_index = jet3.index

            if not local:
                # only computed in skim!
                tree.jet3_BCHMedium = jet3.BCHMedium
                tree.jet3_BCHTight = jet3.BCHTight

            # eta centrality of 3rd leading jet
            tree.jet3_centrality = eventshapes.eta_centrality(
                jet3.fourvect.Eta(),
                jet1.fourvect.Eta(),
                jet2.fourvect.Eta())
            tree.jet3_centrality_boosted = eventshapes.eta_centrality(
                jet3.fourvect_boosted.Eta(),
                jet1.fourvect_boosted.Eta(),
                jet2.fourvect_boosted.Eta())
        elif local:
            # zero the fourvect
            # avoid ghost values from skim if jet selection changed
            tree['jet3_pt'].reset()
            tree['jet3_p'].reset()
            tree['jet3_et'].reset()
            tree['jet3_e'].reset()
            tree['jet3_m'].reset()
            tree['jet3_phi'].reset()
            tree['jet3_eta'].reset()


class TrueTauBlock((TrueTau + MatchedObject).prefix('truetau1_') +
                   (TrueTau + MatchedObject).prefix('truetau2_')):
    dR_truetaus = FloatCol(default=-1)
    dEta_truetaus = FloatCol(default=-1)
    dPhi_truetaus = FloatCol(default=-1)
    theta_truetaus = FloatCol(default=-1)
    cos_theta_truetaus = FloatCol(default=-10)
    truetau_pt_ratio = FloatCol(default=-1)

    @classmethod
    def set(cls, tree, tau1, tau2):
        if tau1.matched:
            truetau = tau1.matched_object
            tree_object = tree.truetau1
            tree_object.nProng = truetau.nProng
            tree_object.nPi0 = truetau.nPi0
            tree_object.charge = truetau.charge
            TrueTau.set(tree_object, truetau.fourvect)
            TrueTau.set_vis(tree_object, truetau.fourvect_vis)
        if tau2.matched:
            truetau = tau2.matched_object
            tree_object = tree.truetau2
            tree_object.nProng = truetau.nProng
            tree_object.nPi0 = truetau.nPi0
            tree_object.charge = truetau.charge
            TrueTau.set(tree_object, truetau.fourvect)
            TrueTau.set_vis(tree_object, truetau.fourvect_vis)
        # angular variables
        if tau1.matched and tau2.matched:
            truetau1 = tau1.matched_object.fourvect_vis
            truetau2 = tau2.matched_object.fourvect_vis
            tree.theta_truetaus = abs(truetau1.Angle(truetau2))
            tree.cos_theta_truetaus = math.cos(tree.theta_truetaus)
            tree.dR_truetaus = truetau1.DeltaR(truetau2)
            tree.dEta_truetaus = abs(truetau2.Eta() - truetau1.Eta())
            tree.dPhi_truetaus = abs(truetau1.DeltaPhi(truetau2))
            # leading pt over subleading pt
            if truetau2.Pt()!=0:
                tree.truetau_pt_ratio = truetau1.Pt() / truetau2.Pt()
            else:
                tree.truetau_pt_ratio = 0

class EventModel(TreeModel):
    trigger = BoolCol(default=True)

    RunNumber = IntCol()
    lbn = IntCol()
    number_of_good_vertices = IntCol()
    averageIntPerXing = FloatCol()
    actualIntPerXing = FloatCol()

    nvtxsoftmet = IntCol()
    nvtxjets = IntCol()

    # event weight given by the PileupReweighting tool
    pileup_weight = FloatCol(default=1.)
    pileup_weight_high = FloatCol(default=1.)
    pileup_weight_low = FloatCol(default=1.)

    mc_weight = FloatCol(default=1.)
    mcevent_pdf_x1_0    = FloatCol(default=1.)
    mcevent_pdf_x2_0    = FloatCol(default=1.)
    mcevent_pdf_id1_0   = FloatCol(default=1.)
    mcevent_pdf_id2_0   = FloatCol(default=1.)
    mcevent_pdf_scale_0 = FloatCol(default=1.)

    #sphericity = FloatCol(default=-1)
    #aplanarity = FloatCol(default=-1)

    #sphericity_boosted = FloatCol(default=-1)
    #aplanarity_boosted = FloatCol(default=-1)

    #sphericity_full = FloatCol(default=-1)
    #aplanarity_full = FloatCol(default=-1)

    sum_pt = FloatCol()
    sum_pt_full = FloatCol()
    vector_sum_pt = FloatCol()
    vector_sum_pt_full = FloatCol()
    resonance_pt = FloatCol()
    true_resonance_pt = FloatCol()
    num_true_jets_no_overlap = IntCol()
    true_jet1_no_overlap_pt = FloatCol(default=-1)
    true_jet2_no_overlap_pt = FloatCol(default=-1)
    true_dEta_jet1_jet2_no_overlap = FloatCol(default=-1)
    true_mass_jet1_jet2_no_overlap = FloatCol(default=-1)
    true_dphi_jj_higgs_no_overlap = FloatCol(default=-9999)

    ntrack_pv = IntCol()
    ntrack_nontau_pv = IntCol()

    # used by the JetCopy filter:
    jet_E_original = stl.vector('float')
    jet_m_original = stl.vector('float')
    jet_pt_original = stl.vector('float')
    jet_eta_original = stl.vector('float')
    jet_phi_original = stl.vector('float')

    error = BoolCol()
    
    @classmethod
    def set(cls, tree, ei):
        tree.RunNumber = ei.runNumber()
        tree.lbn = ei.lumiBlock()

class InclusiveHiggsModel(TreeModel):
    higgs_decay_channel = IntCol(default=-1)


def get_model(datatype, name, prefix=None, is_inclusive_signal=False):
    model = EventModel + MassModel + METModel + RecoTauBlock + RecoJetBlock
    if datatype in (datasets.EMBED, datasets.MCEMBED):
        model += EmbeddingModel
    if datatype != datasets.DATA:
        model += TrueTauBlock
    #if datatype == datasets.MC and 'VBF' in name:
    #    # add branches for VBF Higgs associated partons
    #    model += PartonBlock
    if is_inclusive_signal:
        model += InclusiveHiggsModel
    if prefix is not None:
        return model.prefix(prefix)
    return model

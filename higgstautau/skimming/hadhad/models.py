import ROOT
from rootpy.math.physics.vector import LorentzVector
from rootpy.tree import TreeModel
from rootpy.types import *


class SkimModel(TreeModel):

    number_of_good_vertices = IntCol()
    tau_selected = ROOT.vector('bool')
    tau_numTrack_recounted = ROOT.vector('int')
    pileup_weight = FloatCol(default=1.)
    period_weight = FloatCol(default=1.)
    ggf_weight = FloatCol(default=1.)

    @classmethod
    def reset(cls, tree):

        tree.tau_numTrack_recounted.clear()

    @classmethod
    def set(cls, tree, tau):

        tree.tau_numTrack_recounted.push_back(
            tau.numTrack_recounted)


class TriggerMatching(TreeModel):

    tau_trigger_match_index = ROOT.vector('int')
    tau_trigger_match_thresh = ROOT.vector('int')

    @classmethod
    def reset(cls, tree):

        tree.tau_trigger_match_index.clear()
        tree.tau_trigger_match_thresh.clear()

    @classmethod
    def set(cls, tree, tau):

        tree.tau_trigger_match_index.push_back(
                tau.trigger_match_index)
        tree.tau_trigger_match_thresh.push_back(
                tau.trigger_match_thresh)


class MassModel(TreeModel):

    tau_MMC_mass = FloatCol()
    tau_MMC_resonance = LorentzVector
    tau_MMC_resonance_pt = FloatCol()
    tau_MMC_MET = FloatCol()
    tau_MMC_MET_x = FloatCol()
    tau_MMC_MET_y = FloatCol()
    tau_MMC_MET_phi = FloatCol()

    tau_collinear_mass = FloatCol()
    tau_collinear_momentum_fraction = ROOT.vector('float')

    tau_visible_mass = FloatCol()


class ScaleFactors(TreeModel):

    # tau id efficiency scale factors
    id_eff_sf_loose = ROOT.vector('float')
    id_eff_sf_loose_high = ROOT.vector('float')
    id_eff_sf_loose_low = ROOT.vector('float')

    id_eff_sf_medium = ROOT.vector('float')
    id_eff_sf_medium_high = ROOT.vector('float')
    id_eff_sf_medium_low = ROOT.vector('float')

    id_eff_sf_tight = ROOT.vector('float')
    id_eff_sf_tight_high = ROOT.vector('float')
    id_eff_sf_tight_low = ROOT.vector('float')

    # fakerate scale factors
    fakerate_sf_loose = ROOT.vector('float')
    fakerate_sf_loose_high = ROOT.vector('float')
    fakerate_sf_loose_low = ROOT.vector('float')

    fakerate_sf_medium = ROOT.vector('float')
    fakerate_sf_medium_high = ROOT.vector('float')
    fakerate_sf_medium_low = ROOT.vector('float')

    fakerate_sf_tight = ROOT.vector('float')
    fakerate_sf_tight_high = ROOT.vector('float')
    fakerate_sf_tight_low = ROOT.vector('float')

    # trigger efficiency scale factors
    trigger_eff_sf_loose = ROOT.vector('float')
    trigger_eff_sf_loose_high = ROOT.vector('float')
    trigger_eff_sf_loose_low = ROOT.vector('float')

    trigger_eff_sf_medium = ROOT.vector('float')
    trigger_eff_sf_medium_high = ROOT.vector('float')
    trigger_eff_sf_medium_low = ROOT.vector('float')

    trigger_eff_sf_tight = ROOT.vector('float')
    trigger_eff_sf_tight_high = ROOT.vector('float')
    trigger_eff_sf_tight_low = ROOT.vector('float')


class TauCorrections(ScaleFactors.prefix('tau_')):

    @classmethod
    def reset(cls, tree):

        attrs = ScaleFactors.get_attrs()
        for name, value in attrs:
            getattr(tree.tau, name).clear()

    @classmethod
    def set(cls, tree, tau):

        attrs = ScaleFactors.get_attrs()
        for name, value in attrs:
            getattr(tree.tau, name).push_back(getattr(tau, name))


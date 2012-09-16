import ROOT
from rootpy.math.physics.vector import LorentzVector
from rootpy.tree import TreeModel
from rootpy.types import *


class SkimModel(TreeModel):

    number_of_good_vertices = IntCol()
    tau_selected = ROOT.vector('bool')
    pileup_weight = FloatCol(default=1.)
    ggf_weight = FloatCol(default=1.)


class TriggerMatching(TreeModel):

    tau_trigger_match_index = ROOT.vector('int')
    tau_trigger_match_thresh = ROOT.vector('int')


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


class TauCorrections(TreeModel):

    tau_efficiency_scale_factor = ROOT.vector('float')
    tau_efficiency_scale_factor_high = ROOT.vector('float')
    tau_efficiency_scale_factor_low = ROOT.vector('float')

    tau_fakerate_scale_factor = ROOT.vector('float')
    tau_fakerate_scale_factor_high = ROOT.vector('float')
    tau_fakerate_scale_factor_low = ROOT.vector('float')

    tau_trigger_scale_factor = ROOT.vector('float')
    tau_trigger_scale_factor_high = ROOT.vector('float')
    tau_trigger_scale_factor_low = ROOT.vector('float')

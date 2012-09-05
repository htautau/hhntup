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

    MMC_mass = FloatCol()
    MMC_resonance = LorentzVector
    MMC_resonance_pt = FloatCol()
    MMC_MET = FloatCol()
    MMC_MET_x = FloatCol()
    MMC_MET_y = FloatCol()
    MMC_MET_phi = FloatCol()

    tau_collinear_mass = FloatCol()
    tau_collinear_momentum_fraction = ROOT.vector('float')

    tau_visible_mass = FloatCol()

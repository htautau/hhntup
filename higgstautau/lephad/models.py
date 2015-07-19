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

def get_model(datatype, name, prefix=None, is_inclusive_signal=False):
    model = EventModel
    # if datatype in (datasets.EMBED, datasets.MCEMBED):
    #     model += EmbeddingModel
    # if datatype != datasets.DATA:
    #     model += TrueTauBlock
    #if datatype == datasets.MC and 'VBF' in name:
    #    # add branches for VBF Higgs associated partons
    #    model += PartonBlock
    if prefix is not None:
        return model.prefix(prefix)
    return model

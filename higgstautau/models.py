"""
TreeModels common to both lephad and hadhad analyses
"""

from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy.types import *


class TrueTau_MCBlock(TreeModel):

    visible = BoolCol(default=False)
    hadronic = BoolCol(default=False)
    electron = BoolCol(default=False)
    muon = BoolCol(default=False)
    nprong = IntCol(default=-1111)
    npi0 = IntCol(default=-1111)
    nneutrals = IntCol(default=-1111)
    charge = IntCol()

    fourvect = LorentzVector
    fourvect_vis = LorentzVector
    fourvect_miss = LorentzVector

    fourvect_boosted = LorentzVector
    fourvect_vis_boosted = LorentzVector
    fourvect_miss_boosted = LorentzVector

    prod_vertex = Vector3
    decay_vertex = Vector3
    decay_length = FloatCol(default=-1111)

    dR_tau_nu = FloatCol(default=-1111)
    dTheta3d_tau_nu = FloatCol(default=-1111)


class TrueTau(TreeModel):

    nProng = IntCol(default=-1111)
    nPi0 = IntCol(default=-1111)
    charge = IntCol()

    fourvect = LorentzVector
    fourvect_vis = LorentzVector

    fourvect_boosted = LorentzVector
    fourvect_vis_boosted = LorentzVector


class MatchedObject(TreeModel):

    matched = IntCol()
    matched_dR = FloatCol(default=1111)
    matched_collision = BoolCol()
    matched_pdgId = IntCol()


class Parton(TreeModel):

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector
    pdgId = IntCol()


class PartonBlock((Parton + MatchedObject).prefix('parton1_') + (Parton + MatchedObject).prefix('parton2_')):

    @classmethod
    def set(cls, tree, parton1, parton2):

        # sort by eta
        parton1, parton2 = sorted([parton1, parton2], key=lambda parton: parton.fourvect.Eta())

        parton1_fourvect = parton1.fourvect
        parton2_fourvect = parton2.fourvect

        tree.mass_true_quark1_quark2 = (parton1_fourvect + parton2_fourvect).M()

        tree.parton1_fourvect.set_from(parton1_fourvect)
        tree.parton2_fourvect.set_from(parton2_fourvect)

        beta = (parton1_fourvect + parton2_fourvect).BoostVector()
        tree.parton_beta.set_from(beta)

        parton1_fourvect.Boost(beta * -1)
        parton2_fourvect.Boost(beta * -1)

        tree.parton1_fourvect_boosted.set_from(parton1_fourvect)
        tree.parton2_fourvect_boosted.set_from(parton2_fourvect)

        tree.parton1_pdgId = parton1.pdgId
        tree.parton2_pdgId = parton2.pdgId

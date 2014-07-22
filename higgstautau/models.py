"""
TreeModels common to both lephad and hadhad analyses
"""
from ROOT import TLorentzVector

from rootpy.tree import TreeModel, FloatCol, IntCol, BoolCol
from rootpy.vector import LorentzRotation, LorentzVector, Vector3, Vector2
from rootpy import log
ignore_warning = log['/ROOT.TVector3.PseudoRapidity'].ignore(
    '.*transvers momentum.*')


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


class MMCOutput(FourMomentum.prefix('resonance_')):
    mass = FloatCol()
    MET_et = FloatCol()
    MET_etx = FloatCol()
    MET_ety = FloatCol()
    MET_phi = FloatCol()


class MMCModel(MMCOutput.prefix('mmc0_'),
               MMCOutput.prefix('mmc1_'),
               MMCOutput.prefix('mmc2_')):
    pass


class TrueTau(FourMomentum + FourMomentum.prefix('vis_')):
    nProng = IntCol(default=-1111)
    nPi0 = IntCol(default=-1111)
    charge = IntCol()

    @classmethod
    def set_vis(cls, this, other):
        if isinstance(other, TLorentzVector):
            vect = other
        else:
            vect = other.fourvect
        this.vis_pt = vect.Pt()
        this.vis_p = vect.P()
        this.vis_et = vect.Et()
        this.vis_e = vect.E()
        this.vis_m = vect.M()
        with ignore_warning:
            this.vis_phi = vect.Phi()
            this.vis_eta = vect.Eta()


class MatchedObject(TreeModel):
    matched = BoolCol()
    matched_dR = FloatCol(default=1111)
    #matched_collision = BoolCol()
    #matched_pdgId = IntCol()


class Parton(FourMomentum):
    pdgId = IntCol()


class PartonBlock((Parton + MatchedObject).prefix('parton1_') +
                  (Parton + MatchedObject).prefix('parton2_')):
    dEta_partons = FloatCol(default=-1)
    dR_partons = FloatCol()
    dR_parton_tau = FloatCol()
    parton_beta = Vector3

    @classmethod
    def set(cls, tree, parton1, parton2):
        # TODO: UPDATE

        # sort by eta
        parton1, parton2 = sorted(
            [parton1, parton2],
            key=lambda parton: parton.fourvect.Eta())

        parton1_fourvect = parton1.fourvect
        parton2_fourvect = parton2.fourvect

        tree.mass_true_quark1_quark2 = (
            parton1_fourvect + parton2_fourvect).M()

        tree.parton1_fourvect.copy_from(parton1_fourvect)
        tree.parton2_fourvect.copy_from(parton2_fourvect)

        beta = (parton1_fourvect + parton2_fourvect).BoostVector()
        tree.parton_beta.copy_from(beta)

        parton1_fourvect.Boost(beta * -1)
        parton2_fourvect.Boost(beta * -1)

        tree.parton1_fourvect_boosted.copy_from(parton1_fourvect)
        tree.parton2_fourvect_boosted.copy_from(parton2_fourvect)

        tree.parton1_pdgId = parton1.pdgId
        tree.parton2_pdgId = parton2.pdgId

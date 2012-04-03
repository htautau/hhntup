from rootpy.tree import TreeModel
from rootpy.vector import LorentzVector
from rootpy.types import *
import ROOT

class TrueTau_MCBlock(TreeModel):

    hadronic = BoolCol(default=False)
    nprong = IntCol(default=-1111)
    npi0 = IntCol(default=-1111)
    nneutrals = IntCol(default=-1111)
    charge = IntCol()

    fourvect = LorentzVector
    fourvect_vis = LorentzVector
    fourvect_miss = LorentzVector
    
    fourvect_boost = LorentzVector
    fourvect_vis_boost = LorentzVector
    fourvect_miss_boost = LorentzVector

    dR_tau_nu = FloatCol(default=-1111)
    dTheta3d_tau_nu = FloatCol(default=-1111)


class RecoTau(TreeModel):
    
    BDTJetScore = FloatCol()
    BDTEleScore = FloatCol()
    nPi0 = IntCol()
    seedCalo_numTrack = IntCol()
    charge = IntCol()

    fourvect = LorentzVector
    fourvect_boost = LorentzVector


class TrueTau(TreeModel):
    
    nprong = IntCol(default=-1111)
    npi0 = IntCol(default=-1111)
    charge = IntCol()
    
    fourvect = LorentzVector
    fourvect_vis = LorentzVector


class MatchedObject(TreeModel):

    matched = IntCol()
    matched_dR = FloatCol(default=1111)
    matched_collision = BoolCol()
    matched_pdgId = IntCol()


class EventVariables(TreeModel):

    Mvis_tau1_tau2 = FloatCol()
    numJets = IntCol()
    jet_AntiKt4TopoEM_matched = ROOT.vector("int")
    jet_AntiKt4TopoEM_matched_dR = ROOT.vector("float")
    numVertices = IntCol()
    MET = FloatCol()
    MET_phi = FloatCol()
    HT = FloatCol()
    MMC_mass = FloatCol()
    error = BoolCol()
    mu = IntCol()
    selected = BoolCol()


class RecoJet(TreeModel):
    
    fourvect = LorentzVector
    fourvect_boost = LorentzVector

    jvtxf = FloatCol()


class Parton(TreeModel):

    fourvect = LorentzVector
    pdgId = IntCol()


class RecoTauBlock((RecoTau + MatchedObject).prefix('tau1_') + (RecoTau + MatchedObject).prefix('tau2_')):
    
    @classmethod
    def set(cls, tree, tau1, tau2): pass

 
class RecoJetBlock((RecoJet + MatchedObject).prefix('jet1_') + (RecoJet + MatchedObject).prefix('jet2_')):

    @classmethod
    def set(cls, tree, jet1, jet2): pass


class TrueTauBlock((TrueTau + MatchedObject).prefix('trueTau1_') + (TrueTau + MatchedObject).prefix('trueTau2_')):

    @classmethod
    def set(cls, tree, index, tau): pass


class PartonBlock((Parton + MatchedObject).prefix('parton1_') + (Parton + MatchedObject).prefix('parton2_')):

    @classmethod
    def set(cls, tree, parton1, parton2): pass


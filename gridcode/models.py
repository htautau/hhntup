from rootpy.tree import TreeModel
from rootpy.vector import LorentzVector
from rootpy.types import *
from atlastools.utils import et2pt
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
    
    nProng = IntCol(default=-1111)
    nPi0 = IntCol(default=-1111)
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
    def set(cls, tree, tau1, tau2):

        for i, tau in zip((1,2), (tau1, tau2)):
        
            setattr(tree, 'tau%i_BDTJetScore' % i, tau.BDTJetScore)
            setattr(tree, 'tau%i_BDTEleScore' % i, tau.BDTEleScore)
            setattr(tree, 'tau%i_nPi0' % i, tau.nPi0)
            setattr(tree, 'tau%i_seedCalo_numTrack' % i, tau.seedCalo_numTrack)
            setattr(tree, 'tau%i_charge' % i, tau.charge)
            LorentzVector.__init__(getattr(tree, 'tau%i_fourvect' % i), tau.fourvect)

 
class RecoJetBlock((RecoJet + MatchedObject).prefix('jet1_') + (RecoJet + MatchedObject).prefix('jet2_')):

    @classmethod
    def set(cls, tree, jet1, jet2):
        
        # call copy constructors
        LorentzVector.__init__(tree.jet1_fourvect, jet1.fourvect)
        LorentzVector.__init__(tree.jet2_fourvect, jet2.fourvect)

        try:
            tree.jet1_jvtxf = jet1.jvtxf
            tree.jet2_jvtxf = jet2.jvtxf
        except AttributeError: pass


class TrueTauBlock((TrueTau + MatchedObject).prefix('trueTau1_') + (TrueTau + MatchedObject).prefix('trueTau2_')):

    @classmethod
    def set(cls, tree, index, tau):

        setattr(tree, 'trueTau%i_nProng' % index, tau.nProng)
        setattr(tree, 'trueTau%i_nPi0' % index, tau.nPi0)
        setattr(tree, 'trueTau%i_charge' % index, tau.charge)
    
        getattr(tree, 'trueTau%i_fourvect' % index).SetPtEtaPhiM(
            tau.pt,
            tau.eta,
            tau.phi,
            tau.m)
        
        getattr(tree, 'trueTau%i_fourvect_vis' % index).SetPtEtaPhiM(
            et2pt(tau.vis_Et, tau.vis_eta, tau.vis_m),
            tau.vis_eta,
            tau.vis_phi,
            tau.vis_m)


class PartonBlock((Parton + MatchedObject).prefix('parton1_') + (Parton + MatchedObject).prefix('parton2_')):

    @classmethod
    def set(cls, tree, parton1, parton2):

        # call copy constructors
        LorentzVector.__init__(tree.parton1_fourvect, parton1.fourvect())
        LorentzVector.__init__(tree.parton2_fourvect, parton2.fourvect())

        tree.parton1_pdgId = parton1.pdgId
        tree.parton2_pdgId = parton2.pdgId

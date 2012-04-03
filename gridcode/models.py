from rootpy.tree import TreeModel
from rootpy.vector import LorentzVector, Vector3
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
    
    fourvect_boosted = LorentzVector
    fourvect_vis_boosted = LorentzVector
    fourvect_miss_boosted = LorentzVector

    dR_tau_nu = FloatCol(default=-1111)
    dTheta3d_tau_nu = FloatCol(default=-1111)


class RecoTau(TreeModel):
    
    BDTJetScore = FloatCol()
    BDTEleScore = FloatCol()
    JetBDTLoose = BoolCol()
    JetBDTMedium = BoolCol()
    JetBDTTight = BoolCol()

    nPi0 = IntCol()
    seedCalo_numTrack = IntCol()
    charge = IntCol()

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector


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


class EventVariables(TreeModel):

    Mvis_tau1_tau2 = FloatCol()
    numJets = IntCol()
    jet_fourvect = ROOT.vector('TLorentzVector')
    jet_jvtxf = ROOT.vector('float')
    numVertices = IntCol()
    MET = FloatCol()
    MET_phi = FloatCol()
    HT = FloatCol()
    MMC_mass = FloatCol()
    error = BoolCol()
    mu = IntCol()
    jet_beta = Vector3
    parton_beta = Vector3
    cutflow = IntCol()


class RecoJet(TreeModel):
    
    fourvect = LorentzVector
    fourvect_boosted = LorentzVector

    jvtxf = FloatCol()


class Parton(TreeModel):

    fourvect = LorentzVector
    fourvect_boosted = LorentzVector
    pdgId = IntCol()


class RecoTauBlock((RecoTau + MatchedObject).prefix('tau1_') + (RecoTau + MatchedObject).prefix('tau2_')):
    
    @classmethod
    def set(cls, tree, tau1, tau2):

        for i, tau in zip((1,2), (tau1, tau2)):
        
            fourvect = tau.fourvect
            setattr(tree, 'tau%i_BDTJetScore' % i, tau.BDTJetScore)
            setattr(tree, 'tau%i_BDTEleScore' % i, tau.BDTEleScore)
            
            setattr(tree, 'tau%i_JetBDTLoose' % i, tau.JetBDTLoose)
            setattr(tree, 'tau%i_JetBDTMedium' % i, tau.JetBDTMedium)
            setattr(tree, 'tau%i_JetBDTTight' % i, tau.JetBDTTight)
            
            setattr(tree, 'tau%i_nPi0' % i, tau.nPi0)
            setattr(tree, 'tau%i_seedCalo_numTrack' % i, tau.seedCalo_numTrack)
            setattr(tree, 'tau%i_charge' % i, tau.charge)
            getattr(tree, 'tau%i_fourvect' % i).set_from(fourvect)
            fourvect.Boost(tree.jet_beta * -1)
            getattr(tree, 'tau%i_fourvect_boosted' % i).set_from(fourvect)

 
class RecoJetBlock((RecoJet + MatchedObject).prefix('jet1_') + (RecoJet + MatchedObject).prefix('jet2_')):

    @classmethod
    def set(cls, tree, jet1, jet2):
        
        jet1_fourvect = jet1.fourvect
        jet2_fourvect = jet2.fourvect
        
        tree.jet1_fourvect.set_from(jet1_fourvect)
        tree.jet2_fourvect.set_from(jet2_fourvect)

        beta = (jet1_fourvect + jet2_fourvect).BoostVector()
        tree.jet_beta.set_from(beta)

        tree.jet1_fourvect_boosted.set_from(jet1_fourvect)
        tree.jet2_fourvect_boosted.set_from(jet2_fourvect)

        tree.jet1_fourvect_boosted.Boost(beta * -1)
        tree.jet2_fourvect_boosted.Boost(beta * -1)

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


class PartonBlock((Parton + MatchedObject).prefix('parton1_') + (Parton + MatchedObject).prefix('parton2_')):

    @classmethod
    def set(cls, tree, parton1, parton2):
        
        parton1_fourvect = parton1.fourvect()
        parton2_fourvect = parton2.fourvect()

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

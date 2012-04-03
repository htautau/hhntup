from decorators import memoize
from rootpy.vector import LorentzVector
from rootpy.hep import pdg

"""
This module contains "mixin" classes for adding
functionality to Tree objects ("decorating" them).
"""


class FourMomentum(object):

    def __init__(self):

        self.fourvect_boosted = LorentzVector()
    
    @property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.m)
        return vect


class TauFourMomentum(FourMomentum):

    @property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.seedCalo_eta, self.seedCalo_phi, self.m)
        return vect


class MCParticle(FourMomentum):

    def __init__(self):

        self._particle = pdg.GetParticle(self.pdgId)
        FourMomentum.__init__(self)
    
    def ichildren(self):

        for child in self.child_index:
            try:
                yield getattr(self.tree, self.name)[child]
            except:
                continue

    def traverse_children(self):

        for child in self.ichildren():
            yield child
            for desc in child.traverse_children():
                yield desc
    
    def iparents(self):

        for parent in self.parent_index:
            try:
                yield getattr(self.tree, self.name)[parent]
            except:
                continue
    
    def traverse_parents(self):

        for parent in self.iparents():
            yield parent
            for ancestor in parent.traverse_parents():
                yield ancestor

    def is_leaf(self):

        return not len(self.child_index)

    def final_state(self):

        if self.is_leaf():
            return [self]
        return [particle for particle in self.traverse_children() if particle.is_leaf()]
    
    def fourvect(self):
        
        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self._particle.Mass()*1000.)
        return vect 
    
    def __repr__(self):

        return self.__str__()
    
    def __str__(self):

        return "%s (m: %.3f MeV, pt: %.1f GeV, eta: %.2f, phi: %.2f)" % \
            (self._particle.GetName(), self._particle.Mass() * 1000., self.pt / 1000., self.eta, self.phi)

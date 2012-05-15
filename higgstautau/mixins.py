from decorators import cached_property
from rootpy.math.physics.vector import LorentzVector
from rootpy.hep import pdg
from atlastools.units import GeV

"""
This module contains "mixin" classes for adding
functionality to Tree objects ("decorating" them).
"""


class FourMomentum(object):

    def __init__(self):

        self.fourvect_boosted = LorentzVector()

    @cached_property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.m)
        return vect

    def __repr__(self):

        return self.__str__()

    def __str__(self):

        return "%s (m: %.3f MeV, pt: %.1f MeV, eta: %.2f, phi: %.2f)" % \
            (self.__class__.__name__,
             self.m,
             self.pt,
             self.eta,
             self.phi)


class TauFourMomentum(FourMomentum):

    def __init__(self):

        self.weight = 1.

        self.centrality = 0.
        self.centrality_boosted = 0.

        self.matched = False
        self.matched_dR = 9999.
        self.matched_collision = False

        super(TauFourMomentum, self).__init__()

    @cached_property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.m)
        return vect

    @cached_property
    def leadtrack_idx(self):
        """
        Return index of leading track
        """
        ldtrkindex = -1
        ldtrkpt = 0.
        for i in xrange(self.track_n):
            pt = self.track_pt[i]
            if pt > ldtrkpt:
                ldtrkindex = i
                ldtrkpt = pt
        return ldtrkindex



class ElFourMomentum(FourMomentum):

    @cached_property
    def fourvect(self):

        e   = self.cl_E
        eta = self.tracketa
        phi = self.trackphi
        et  = e/cosh(eta)

        vect = LorentzVector()
        vect.SetPtEtaPhiM(et, eta, phi, self.m)
        return vect


class MCTauFourMomentum(FourMomentum):

    @cached_property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self.m)
        return vect


class MCParticle(FourMomentum):

    def __init__(self):

        self._particle = pdg.GetParticle(self.pdgId)
        FourMomentum.__init__(self)

    def ichildren(self):

        try:
            for child in self.child_index:
                yield getattr(self.tree, self.name)[child]
        except GeneratorExit:
            pass

    def traverse_children(self):

        try:
            for child in self.ichildren():
                yield child
                for desc in child.traverse_children():
                    yield desc
        except GeneratorExit:
            pass

    def iparents(self):

        try:
            for parent in self.parent_index:
                yield getattr(self.tree, self.name)[parent]
        except GeneratorExit:
            pass

    def traverse_parents(self):

        try:
            for parent in self.iparents():
                yield parent
                for ancestor in parent.traverse_parents():
                    yield ancestor
        except GeneratorExit:
            pass

    def is_leaf(self):

        return not len(self.child_index)

    @cached_property
    def first_self(self):

        for parent in self.iparents():
            if parent.pdgId == self.pdgId:
                return parent.first_self
        return self

    @cached_property
    def last_self(self):

        for child in self.ichildren():
            if child.pdgId == self.pdgId:
                return child.last_self
        return self

    @cached_property
    def final_state(self):

        if self.is_leaf():
            return [self]
        return [particle for particle in self.traverse_children() if particle.is_leaf()]

    @cached_property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(self.pt, self.eta, self.phi, self._particle.Mass() * GeV)
        return vect

    def __repr__(self):

        return self.__str__()

    def __str__(self):

        return ("%s ("
                "status: %d, "
                "m: %.3f MeV, pt: %.1f GeV, eta: %.2f, phi: %.2f, "
                "x: %.4f, y: %.4f, z: %.4f)") % \
            (self._particle.GetName(),
             self.status,
             self._particle.Mass() * GeV,
             self.pt / GeV,
             self.eta, self.phi,
             self.vx_x,
             self.vx_y,
             self.vx_z)

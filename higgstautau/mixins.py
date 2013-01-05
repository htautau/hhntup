import math
from decorators import cached_property
from rootpy.math.physics.vector import LorentzVector
from rootpy.extern.hep import pdg
from atlastools.units import GeV

"""
This module contains "mixin" classes for adding
functionality to Tree objects ("decorating" them).
"""

__all__ = [
    'FourMomentum',
    'TauFourMomentum',
    'TauFourMomentumSkim',
    'ElectronFourMomentum',
    'MCTauFourMomentum',
    'MCParticle',
]


class MatchedObject(object):

    def __init__(self):

        self.matched = False
        self.matched_dR = 9999.
        self.matched_collision = False


class FourMomentum(MatchedObject):

    def __init__(self):

        self.fourvect_boosted = LorentzVector()
        super(FourMomentum, self).__init__()

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

        super(TauFourMomentum, self).__init__()

        self.centrality = 0.
        self.centrality_boosted = 0.

        # vertex association
        self.vertex_prob = 0.

        # overlap checking
        self.min_dr_jet = 9999.

        self._pt_nominal = -1111.

    @property
    def pt_nominal(self):

        if self._pt_nominal != -1111.:
            return self._pt_nominal
        return self.pt

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


class TauFourMomentumSkim(TauFourMomentum):

    def __init__(self):

        super(TauFourMomentumSkim, self).__init__()

        SF_DEFAULT = -1.

        self.trigger_match_thresh = 0
        self.trigger_match_index = -1

        # tau id efficiency scale factors
        self.id_eff_sf_loose = SF_DEFAULT
        self.id_eff_sf_loose_high = SF_DEFAULT
        self.id_eff_sf_loose_low = SF_DEFAULT

        self.id_eff_sf_medium = SF_DEFAULT
        self.id_eff_sf_medium_high = SF_DEFAULT
        self.id_eff_sf_medium_low = SF_DEFAULT

        self.id_eff_sf_tight = SF_DEFAULT
        self.id_eff_sf_tight_high = SF_DEFAULT
        self.id_eff_sf_tight_low = SF_DEFAULT

        # fakerate scale factors
        self.fakerate_sf_loose = SF_DEFAULT
        self.fakerate_sf_loose_high = SF_DEFAULT
        self.fakerate_sf_loose_low = SF_DEFAULT

        self.fakerate_sf_medium = SF_DEFAULT
        self.fakerate_sf_medium_high = SF_DEFAULT
        self.fakerate_sf_medium_low = SF_DEFAULT

        self.fakerate_sf_tight = SF_DEFAULT
        self.fakerate_sf_tight_high = SF_DEFAULT
        self.fakerate_sf_tight_low = SF_DEFAULT

        # trigger efficiency scale factors
        self.trigger_eff_sf_loose = SF_DEFAULT
        self.trigger_eff_sf_loose_high = SF_DEFAULT
        self.trigger_eff_sf_loose_low = SF_DEFAULT

        self.trigger_eff_sf_medium = SF_DEFAULT
        self.trigger_eff_sf_medium_high = SF_DEFAULT
        self.trigger_eff_sf_medium_low = SF_DEFAULT

        self.trigger_eff_sf_tight = SF_DEFAULT
        self.trigger_eff_sf_tight_high = SF_DEFAULT
        self.trigger_eff_sf_tight_low = SF_DEFAULT

        # colliniear mass approx
        self.collinear_momentum_fraction = -9999.

        # track recounting
        self.numTrack_recounted = -1


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

    @cached_property
    def num_children(self):

        return len(self.child_index)

    @cached_property
    def num_parents(self):

        return len(self.parent_index)

    def get_child(self, index):

        index = self.child_index[index]
        return getattr(self.tree, self.name)[index]

    def get_parent(self, index):

        index = self.parent_index[index]
        return getattr(self.tree, self.name)[index]

    def iter_children(self):

        try:
            for child in self.child_index:
                yield getattr(self.tree, self.name)[child]
        except GeneratorExit:
            pass

    def iter_parents(self):

        try:
            for parent in self.parent_index:
                yield getattr(self.tree, self.name)[parent]
        except GeneratorExit:
            pass

    def traverse_children(self):

        try:
            for child in self.iter_children():
                yield child
                for desc in child.traverse_children():
                    yield desc
        except GeneratorExit:
            pass

    def traverse_parents(self):

        try:
            for parent in self.iter_parents():
                yield parent
                for ancestor in parent.traverse_parents():
                    yield ancestor
        except GeneratorExit:
            pass

    def is_stable(self):

        return self.status == 1

    @cached_property
    def first_self(self):

        for parent in self.iter_parents():
            if parent.pdgId == self.pdgId:
                return parent.first_self
        return self

    @cached_property
    def last_self(self):

        for child in self.iter_children():
            if child.pdgId == self.pdgId:
                return child.last_self
        return self

    @cached_property
    def final_state(self):

        if self.is_stable():
            return [self]
        return [particle for particle in self.traverse_children()
                if particle.is_stable()]

    @cached_property
    def fourvect(self):

        vect = LorentzVector()
        vect.SetPtEtaPhiM(
                self.pt,
                self.eta,
                self.phi,
                self.m)
        #        self._particle.Mass() * GeV)
        return vect

    def export_graphvis(self, out_file=None):

        def particle_to_str(particle):

            return ('%s\\n'
                    'mass = %.3f MeV\\n'
                    'pt = %.3f GeV\\n'
                    'eta = %.2f\\n'
                    'status = %d') % (
                    particle._particle.GetName(),
                    #particle._particle.Mass() * GeV,
                    particle.m,
                    particle.pt / GeV,
                    particle.eta,
                    particle.status)

        def recurse(particle, parent=None):

            out_file.write('%d [label="%s"] ;\n' % (
                particle.barcode, particle_to_str(particle)))

            if parent is not None:
                # Add edge to parent
                out_file.write('%d -> %d ;\n' % (
                    parent.barcode,
                    particle.barcode))

            # recurse on children
            for child in particle.iter_children():
                recurse(child, particle)

        close_file = True
        if out_file is None:
            out_file = open('event.dot', 'w')
        elif isinstance(out_file, basestring):
            out_file = open(out_file, 'w')
        else:
            close_file = False

        out_file.write('digraph Tree {\n')
        out_file.write('size="7.5,10" ;\n')
        out_file.write('orientation=landscape ;\n')
        recurse(self, None)
        out_file.write('}')

        if close_file:
            out_file.close()
        else:
            return out_file

    def __repr__(self):

        return self.__str__()

    def __str__(self):

        return ("%s ("
                "status: %d, "
                "m: %.3f MeV, pt: %.1f GeV, eta: %.2f, phi: %.2f, "
                "x: %.4f, y: %.4f, z: %.4f)") % \
            (self._particle.GetName(),
             self.status,
             #self._particle.Mass() * GeV,
             self.m,
             self.pt / GeV,
             self.eta, self.phi,
             self.vx_x,
             self.vx_y,
             self.vx_z)


class ElectronFourMomentum(FourMomentum):

    @cached_property
    def fourvect(self):

        if ((self.nSCTHits + self.nPixHits) < 4):
            # electron with low number of tracker hits
            eta = self.cl_eta
            phi = self.cl_phi
            et  = self.cl_E / math.cosh(self.cl_eta)
        else:
            eta = self.tracketa
            phi = self.trackphi
            et  = self.cl_E / math.cosh(self.tracketa)

        vect = LorentzVector()
        vect.SetPtEtaPhiE(et, eta, phi, self.cl_E)
        return vect

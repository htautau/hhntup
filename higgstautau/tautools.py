from rootpy.vector import Vector3, LorentzVector as FourVector
from rootpy.extern.hep import pdg
from rootpy import asrootpy
from .decorators import cached_property

try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

"""
This module contains utilities for working with
MC tau decays in the MC block of xAODs
"""

class TauDecay(object):

    def __init__(self, initial_state):
        # get tau just before decay
        self.init = initial_state
        self.valid = True
        self.vertex = self.init.decayVtx()
        # traverse to the final state while counting unique particles
        # first build list of unique children
        # (ignore copies in the event record)
        children = [
            self.vertex.outgoingParticle(i) for i in 
            xrange(self.vertex.nOutgoingParticles())
            ]
        stable_children = []
        for child in children:
            if child.status() == 1 or child.status() == 2:
                stable_children.append(child)

        # classify final state particles
        neutrinos = []
        charged_pions = []
        charged_kaons = []
        neutral_kaons = []
        electrons = []
        muons = []
        photons = []

        hadronic = False
        leptonic_electron = False
        leptonic_muon = False
        nprong = 0

        child_pdgid_freq = {}
        for p in stable_children:
            pdgid = p.absPdgId()
            if pdgid not in child_pdgid_freq:
                child_pdgid_freq[pdgid] = 1
            else:
                child_pdgid_freq[pdgid] += 1
            if p.isNeutrino():
                neutrinos.append(p)
                if pdgid == pdg.nu_mu:
                    leptonic_muon = True
                elif pdgid == pdg.nu_e:
                    leptonic_electron = True
            if pdgid == pdg.pi_plus:
                hadronic = True
                nprong += 1
                charged_pions.append(p)
            elif pdgid == pdg.gamma:
                photons.append(p)
            elif pdgid == pdg.e_minus:
                electrons.append(p)
            elif pdgid == pdg.mu_minus:
                muons.append(p)
            elif pdgid == pdg.K_plus:
                hadronic = True
                nprong += 1
                charged_kaons.append(p)
            elif pdgid in (pdg.K_S0, pdg.K_L0, pdg.K0):
                neutral_kaons.append(p)

        self.child_pdgid_freq = child_pdgid_freq
        self.neutrinos = neutrinos
        self.charged_pions = charged_pions
        self.charged_kaons = charged_kaons
        self.neutral_kaons = neutral_kaons
        self.electrons = electrons
        self.muons = muons
        self.photons = photons

        self.hadronic = hadronic
        self.leptonic_electron = leptonic_electron
        self.leptonic_muon = leptonic_muon
        self.nprong = nprong

        self.matched = False
        self.matched_object = None


    @property
    def privtx(self):
        return self.vertex

    @property
    def privtx_x(self):
        return self.vertex.x()

    @property
    def privtx_y(self):
        return self.vertex.y()

    @property
    def privtx_z(self):
        return self.vertex.z()

    @cached_property
    def decay_vertex(self):
        nu_tau = None
        nu_tau_vertex = None
        for nu in self.neutrinos:
            if nu.absPdgId() == pdg.nu_tau:
                nu_tau = nu
                nu_tau_vertex = nu.decayVtx()
                break

        if nu_tau is None:
            return Vector3(0, 0, 0)
        
        return Vector3(nu_tau_vertex.x(),
                       nu_tau_vertex.y(),
                       nu_tau_vertex.z())
    @property
    def secvtx(self):
        return self.decay_vertex

    @property
    def secvtx_x(self):
        return self.decay_vertex.X()

    @property
    def secvtx_y(self):
        return self.decay_vertex.Y()

    @property
    def secvtx_z(self):
        return self.decay_vertex.Z()

    @cached_property
    def decay_vect(self):
        return self.decay_vertex - self.prod_vertex

    @cached_property
    def decay_length(self):
        return self.decay_vect.Mag()

    @cached_property
    def decay_angle(self):
        return self.decay_vect.Angle(self.fourvect_visible)

    @cached_property
    def npi0(self):
        if pdg.pi0 in self.child_pdgid_freq:
            return self.child_pdgid_freq[pdg.pi0]
        return 0

    @cached_property
    def nneutrals(self):
        return self.npi0 + len(self.neutral_kaons)

    @cached_property
    def charge(self):
        return self.init.charge()

    @cached_property
    def fourvect(self):
        return asrootpy(self.init.p4())

    @cached_property
    def fourvect_visible(self):
        return self.fourvect - self.fourvect_missing

    @cached_property
    def fourvect_missing(self):
        missing = FourVector()
        return missing + sum([asrootpy(p.p4()) for p in self.neutrinos])

    @cached_property
    def dr_vistau_nu(self):
        return self.fourvect_visible.DeltaR(self.fourvect_missing)

    @cached_property
    def dtheta3d_vistau_nu(self):
        return self.fourvect_visible.Angle(self.fourvect_missing)

    def __str__(self):
        return self.__repr__()

    def __repr__(self):
        output = StringIO.StringIO()
        print >> output, "initial state:"
        print >> output, self.init
        print >> output, 'npi0: %d' % self.npi0
        print >> output, 'nprong: %d' % self.nprong
        # print >> output, "final state:"
        # for thing in self.final:
        #     print >> output, " - %s" % thing
        rep = output.getvalue()
        output.close()
        return rep


def get_particles(event, pdgid, num_expected=None):
    if not isinstance(pdgid, (list, tuple)):
        pdgid = [pdgid]
    particles = []
    for mc in event.mc:
        if mc.pdgId in pdgid:
            particles.append(mc)
            if num_expected is not None and len(particles) == num_expected:
                break
    return particles


def get_tau_decays(truth_taus, num_expected=None):
    """
    Get all taus and their decay products
    """
    decays = []
    for tau in truth_taus:
        decays.append(TauDecay(tau))
        if num_expected is not None and len(decays) == num_expected:
            break
    return decays

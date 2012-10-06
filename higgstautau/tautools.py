from rootpy.math.physics.vector import Vector3, LorentzVector as FourVector
from rootpy.hep import pdg

from .decorators import cached_property

try:
    import cStringIO as StringIO
except ImportError:
    import StringIO

"""
This module contains utilities for working with
MC tau decays in the MC block of D3PDs.
"""

class TauDecay(object):

    def __init__(self, initial_state):

        # get tau just before decay
        self.init = initial_state.last_self
        # traverse to the final state while counting unique particles
        # first build list of unique children
        # (ignore copies in the event record)
        children = list(set([
            child.last_self for child in
            self.init.traverse_children()]))
        # count the frequency of each pdgId
        child_pdgid_freq = {}
        for child in children:
            pdgid = abs(child.pdgId)
            if pdgid not in child_pdgid_freq:
                child_pdgid_freq[pdgid] = 1
            else:
                child_pdgid_freq[pdgid] += 1
        self.children = children
        self.child_pdgid_freq = child_pdgid_freq
        # collect particles in the final state
        self.final = [p for p in children if p.is_leaf()]
        # some decays are not fully stored in the D3PDs
        # flag them...
        self.complete = True
        if len(self.final) <= 1:
            self.complete = False
        # classify final state particles
        neutrinos = []
        charged_pions = []
        charged_kaons = []
        neutral_kaons = []
        electrons = []
        muons = []
        photons = []
        for p in self.final:
            if abs(p.pdgId) in (pdg.nu_e, pdg.nu_mu, pdg.nu_tau):
                neutrinos.append(p)
            if abs(p.pdgId) == pdg.pi_plus:
                charged_pions.append(p)
            elif p.pdgId == pdg.gamma:
                photons.append(p)
            elif abs(p.pdgId) == pdg.e_minus:
                electrons.append(p)
            elif abs(p.pdgId) == pdg.mu_minus:
                muons.append(p)
            elif abs(p.pdgId) == pdg.K_plus:
                charged_kaons.append(p)
            elif p.pdgId in (pdg.K_S0, pdg.K_L0, pdg.K0):
                neutral_kaons.append(p)
        self.neutrinos = neutrinos
        self.charged_pions = charged_pions
        self.charged_kaons = charged_kaons
        self.neutral_kaons = neutral_kaons
        self.electrons = electrons
        self.muons = muons
        self.photons = photons

    @cached_property
    def has_charged_rho(self):

        if pdg.rho_plus in self.child_pdgid_freq:
            return True
        return False

    @cached_property
    def has_a1(self):

        if pdg.a_1_plus in self.child_pdgid_freq:
            return True
        return False

    @cached_property
    def prod_vertex(self):

        return Vector3(self.init.vx_x,
                       self.init.vx_y,
                       self.init.vx_z)

    @cached_property
    def decay_vertex(self):

        nu_tau = None
        # use production vertex of nu_tau
        last_tau = self.init.last_self
        for child in last_tau.iter_children():
            if abs(child.pdgId) == pdg.nu_tau:
                nu_tau = child
                break
        if nu_tau is None:
            return Vector3(0, 0, 0)
        return Vector3(nu_tau.vx_x,
                       nu_tau.vx_y,
                       nu_tau.vx_z)

    @cached_property
    def decay_length(self):

        return (self.decay_vertex - self.prod_vertex).Mag()

    @cached_property
    def npi0(self):

        if pdg.pi0 in self.child_pdgid_freq:
            return self.child_pdgid_freq[pdg.pi0]
        return 0

    @cached_property
    def electron(self):
        """
        Return True if this is a decay to an electron
        """
        return len(self.electrons) > 0

    @cached_property
    def muon(self):
        """
        Return True if this is a decay to a muon
        """
        return len(self.muons) > 0

    @cached_property
    def hadronic(self):
        """
        Return True if this is a hadronic decay else False for leptonic
        """
        return any(self.charged_pions + self.charged_kaons)

    @cached_property
    def nprong(self):
        """
        Return number of charged particles in final state
        (for hadronic decays only)
        """
        return len(self.charged_pions + self.charged_kaons)

    @cached_property
    def nneutrals(self):

        return self.npi0 + len(self.neutral_kaons)

    @cached_property
    def charge(self):

        return self.init.charge

    @cached_property
    def fourvect(self):

        return self.init.fourvect

    @cached_property
    def fourvect_visible(self):

        return self.fourvect - self.fourvect_missing

    @cached_property
    def fourvect_missing(self):

        missing = FourVector()
        return missing + sum([p.fourvect for p in self.neutrinos])

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
        print >> output, "final state:"
        for thing in self.final:
            print >> output, " - %s" % thing
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


def get_tau_decays(event, parent_pdgid=None, status=None, num_expected=None):
    """
    Get all taus and their decay products

    parent_pdgid: pdgid or list of pdgids of accepted parent particles
    status: accepted status
    """
    if parent_pdgid is not None:
        if not isinstance(parent_pdgid, (list, tuple)):
            parent_pdgid = [parent_pdgid]
    if status is not None:
        if not isinstance(status, (list, tuple)):
            status = [status]
    else:
        # 2 for Pythia, 11 for Herwig, 195 for AlpgenJimmy
        status=(2, 11, 195)
    decays = []
    # TODO speed this up by recursing from primary interaction
    for mc in event.mc:
        if mc.pdgId in (pdg.tau_plus, pdg.tau_minus) and mc.status in status:
            if parent_pdgid is not None:
                accept = False
                for parent in mc.first_self.iter_parents():
                    if parent.pdgId in parent_pdgid:
                        accept = True
                        break
                if not accept:
                    continue
            decays.append(TauDecay(mc))
            if num_expected is not None and len(decays) == num_expected:
                break
    return decays

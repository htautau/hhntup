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

    def __init__(self, initial_state, final_state):
        """
        A tau decay (for all intents and purposes here)
        is composed of an initial and final state.
        """
        self.init = initial_state
        self.final = final_state
        # some decays are not fully stored in the D3PDs
        # flag them...
        self.complete = True
        if len(final_state) == 1:
            self.complete = False

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
        for child in last_tau.ichildren():
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

        npi0 = 0
        for child in self.init.traverse_children():
            if child.pdgId == pdg.pi0:
                npi0 += 1
        return npi0

    @cached_property
    def charged_pions(self):
        """
        Return all charged pions in final state
        """
        return [p for p in self.final if p.pdgId in (pdg.pi_minus, pdg.pi_plus)]

    @cached_property
    def charged_kaons(self):

        return [p for p in self.final if p.pdgId in (pdg.K_minus, pdg.K_plus)]

    @cached_property
    def neutral_kaons(self):

        return [p for p in self.final if p.pdgId in (pdg.K_S0, pdg.K_L0)]

    @cached_property
    def photons(self):

        return [p for p in self.final if p.pdgId == pdg.gamma]

    @cached_property
    def neutrinos(self):
        """
        Return all neutrinos in final state
        """
        return [p for p in self.final if abs(p.pdgId) in (pdg.nu_e, pdg.nu_mu, pdg.nu_tau)]

    @cached_property
    def electrons(self):
        """
        Return all electrons in final state
        """
        return [p for p in self.final if abs(p.pdgId) == pdg.e_minus]

    @cached_property
    def muons(self):
        """
        Return all muons in final state
        """
        return [p for p in self.final if abs(p.pdgId) == pdg.mu_minus]

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
    def fourvect(self):

        return self.init.fourvect

    @cached_property
    def fourvect_visible(self):

        #return sum([p.fourvect for p in self.charged_pions + self.photons + self.charged_kaons + self.neutral_kaons])
        return self.fourvect - self.fourvect_missing

    @cached_property
    def fourvect_missing(self):

        missing = FourVector()
        return missing + sum([p.fourvect for p in self.neutrinos])

    @cached_property
    def dR_tau_nu(self):

        return self.fourvect_visible.DeltaR(self.fourvect_missing)

    @cached_property
    def dTheta3d_tau_nu(self):

        return self.fourvect_visible.Angle(self.fourvect_missing)

    def __str__(self):

        return self.__repr__()

    def __repr__(self):

        output = StringIO.StringIO()
        print >> output, self.init
        for thing in self.final:
            print >> output, "\t%s" % thing
        rep = output.getvalue()
        output.close()
        return rep


def get_tau_decays(event, parent_pdgid=None, status=(2, 11)):
    """
    Get all taus and their decay products

    parent_pdgid: pdgid or list of pdgids of accepted parent particles
    status: accepted status (include 11 for Herwig)
    """
    if parent_pdgid is not None:
        if not isinstance(parent_pdgid, (list, tuple)):
            parent_pdgid = [parent_pdgid]
    if status is not None:
        if not isinstance(status, (list, tuple)):
            status = [status]
    decays = []
    for mc in event.mc:
        if mc.pdgId in (pdg.tau_plus, pdg.tau_minus):
            if status is None or mc.status in status:
                init_state = mc
                if parent_pdgid is not None:
                    accept = True
                    for parent in mc.iparents():
                        if parent.pdgId not in parent_pdgid:
                            accept = False
                            break
                    if not accept:
                        continue
                final_state = mc.final_state
                decays.append(TauDecay(init_state, final_state))
    return decays

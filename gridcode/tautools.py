from decorators import memoize
from rootpy.hep.vector import FourVector
from rootpy.hep import pdg

try:
    import cStringIO as StringIO
except ImportError:
    import StringIO


class TauDecay(object):

    def __init__(self, initial_state, final_state):
        """
        A tau decay (for all intents and purposes here)
        is composed of an initial and final state
        """
        self.init = initial_state
        self.final = final_state
        # some decays are not fully stored in the D3PDs
        # flag them...
        self.complete = True
        if len(final_state) == 1:
            self.complete = False
    
    @property
    def charged_pions(self):
        """
        Return all charged pions in final state
        """
        return [p for p in self.final if p.pdgId in (pdg.pi_minus, pdg.pi_plus)]
    
    @property
    def charged_kaons(self):

        return [p for p in self.final if p.pdgId in (pdg.K_minus, pdg.K_plus)]
    
    @property
    def neutral_kaons(self):
        
        return [p for p in self.final if p.pdgId in (pdg.K_S0, pdg.K_L0)]
    
    @property
    def photons(self):

        return [p for p in self.final if p.pdgId == pdg.gamma]

    @property
    def neutrinos(self):
        """
        Return all neutrinos in final state
        """
        return [p for p in self.final if abs(p.pdgId) in (pdg.nu_e, pdg.nu_mu, pdg.nu_tau)]
    
    @property
    def hadronic(self):
        """
        Return True if this is a hadronic decay else False for leptonic
        """
        return any(self.charged_pions + self.charged_kaons)
    
    @property
    def nprong(self):
        """
        Return number of charged particles in final state
        """
        return len(self.charged_pions + self.charged_kaons)
    
    @property
    def npi0(self):
        """
        Return number of neutral pions: #(gamma)/2
        """
        return len(self.photons) / 2

    @property
    def nneutrals(self):

        return self.npi0 + len(self.neutral_kaons)
    
    @property
    def fourvect(self):
        
        return self.init.fourvect()
    
    @property
    def fourvect_visible(self):

        return sum([p.fourvect() for p in self.charged_pions + self.photons + self.charged_kaons + self.neutral_kaons])

    @property
    def fourvect_missing(self):

        return sum([p.fourvect() for p in self.neutrinos])

    @property
    def dR_tau_nu(self):

        return self.fourvect.DeltaR(self.fourvect_missing)
        
    @property
    def dTheta3d_tau_nu(self):

        return 0.

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


def get_tau_decays(event):
    """
    Get all taus and their decay products
    """
    decays = []
    for mc in event.mc:
        if mc.pdgId in (pdg.tau_plus, pdg.tau_minus) and mc.status == 2:
            init_state = mc
            final_state = mc.final_state()
            decays.append(TauDecay(init_state, final_state))
    return decays

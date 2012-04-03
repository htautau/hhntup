from atlastools import pdg
from decorators import memoize
from rootpy.hep.vector import FourVector

class TauDecay(object):

    def __init__(self, initial_state, final_state):
        """
        A tau decay (for all intents and purposes here)
        is composed of an initial and final state
        """
        self.init = initial_state
        self.final = final_state
    
    @property
    @memoize
    def charged_pions(self):
        """
        Return all charged pions in final state
        """
        return [p in self.final if p.pdgId in (pdg.PI_PLUS, pdg.PI_MINUS)]
    
    @property
    @memoize
    def photons(self):

        return [p in self.final if p.pdgId == pdg.GAMMA]

    @property
    @memoize
    def neutrinos(self):
        """
        Return all neutrinos in final state
        """
        return [p in self.final if p.pdgId in (pdg.NUE, pdg.NUM, pdg.NUTAU)]
    
    @property
    @memoize
    def hadronic(self):
        """
        Return True if this is a hadronic decay else False for leptonic
        """
        return any(self.charged_pions)
    
    @property
    @memoize
    def nprong(self):
        """
        Return number of charged particles in final state
        """
        return len(self.charged_pions)
    
    @property
    @memoize
    def npi0(self):
        """
        Return number of neutral pions: #(gamma)/2
        """
        return len(self.photons) / 2
    
    @property
    def fourvect(self):
        
        return self.init.fourvect 
    
    @property
    @memoize
    def fourvect_visible(self):

        return sum([p.fourvect for p in self.charged_pions + self.photons])

    @property
    @memoize
    def fourvect_missing(self):

        return sum([p.fourvect for p in self.neutrinos])

    @property
    @memoize
    def dR_tau_nu(self):

        return self.fourvect.DeltaR(self.fourfect_missing)
        
    @property
    @memoize
    def dTheta3d_tau_nu(self):

        return 0.

    @memoize
    def __str__(self):

        return self.__repr__()
    
    @memoize
    def __repr__(self):

        print self.init
        for thing in self.final:
            print "\t%s" % thing


def get_tau_decays(event):
    """
    Get all taus and their decay products
    """
    decays = []
    for mc in event.mc:
        if mc.pdgId in (pdg.tau_plus, pdg.tau_minus):
            init_state = mc
            final_state = mc.final_state()
            # some decays are not fully stored in the D3PDs
            # ignore them...
            if len(final_state) > 1:
                decays.append(TauDecay(init_state, final_state))
    return decays

def get_VBF_partons(event):
    
    return []

from atlastools import pdg

class MCParticle(object):

    def __init__(self):

        self.__particle = pdg.GetParticle(self.pdgId)
    
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

    def __repr__(self):

        return self.__str__()
    
    def __str__(self):

        return "%s (m: %.3f MeV, pt: %.1f GeV, eta: %.2f, phi: %.2f)" % (self.__particle.GetName(), self.__particle.Mass()*1000., self.pt/1000., self.eta, self.phi)

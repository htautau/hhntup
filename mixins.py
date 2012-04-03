
class MCParticle(object):

    def iterchildren(self):

        for child in self.children(): 
            yield getattr(self.tree, self.name)[child]

    def iterparents(self):

        for parent in self.parents(): 
            yield getattr(self.tree, self.name)[parent]


class MCParticle(object):

    def iterchildren(self):

        for child in self.children: 
            try:
                yield getattr(self.tree, self.name)[child]
            except:
                continue

    def iterparents(self):

        for parent in self.parents: 
            try:
                yield getattr(self.tree, self.name)[parent]
            except:
                continue

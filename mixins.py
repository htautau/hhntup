
class MCParticle(object):

    def ichildren(self):

        for child in self.children:
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

        for parent in self.parents:
            try:
                yield getattr(self.tree, self.name)[parent]
            except:
                continue
    
    def traverse_parents(self):

        for parent in self.iparents():
            yield parent
            for ancestor in parent.traverse_parents():
                yield ancestor

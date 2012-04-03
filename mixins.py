
class MCParticle(object):

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

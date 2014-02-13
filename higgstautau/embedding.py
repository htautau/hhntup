from rootpy.tree.filtering import EventFilter
from . import log; log = log[__name__]


class EmbeddingPileupPatch(EventFilter):

    def passes(self, event):
        # fix averageIntPerXing
        # HACK: stored in pz of mc particle with pdg ID 39
        # https://twiki.cern.ch/twiki/bin/viewauth/Atlas/EmbeddingTools
        averageIntPerXing = None
        for p in event.mc:
            if p.pdgId == 39:
                averageIntPerXing = p.fourvect.Pz()
                break
        if averageIntPerXing is not None:
            event.averageIntPerXing = averageIntPerXing
        else:
            log.warning("pdgID 39 not found! Skipping event...")
            # ignore event
            return None
        return True


class EmbeddingIsolation(EventFilter):
    """
    2012 embedding isolation systematics
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/EmbeddingTools
    """
    def __init__(self, tree, **kwargs):
        self.tree = tree
        super(EmbeddingIsolation, self).__init__(**kwargs)

    def passes(self, event):
        # no isolation
        isolation = 0
        found = False
        # find particle with pdgid = 82
        for p in event.mc:
            if p.pdgId == 82:
                found = True
                if p.fourvect.Px() > 0:
                    # default isolation
                    isolation = 1
                    if p.fourvect.Py() > 0:
                        # tight isolation
                        isolation = 2
                break
        if not found:
            log.warning("pdgID 82 not found! Skipping event...")
            # ignore event
            return None
        self.tree.embedding_isolation = isolation
        return True

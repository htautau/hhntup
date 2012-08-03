from rootpy.tree.filtering import EventFilter
from atlastools import utils


class TauJetOverlapRemoval(EventFilter):
    """
    Precendence: taus > jets
    """
    def __init__(self, dr=.2, **kwargs):

        super(TauJetOverlapRemoval, self).__init__(**kwargs)
        self.dr = dr

    def passes(self, event):

        # remove overlap with taus
        event.jets.select(lambda jet:
                not any([tau for tau in event.taus if
                (utils.dR(jet.eta, jet.phi, tau.eta, tau.phi) < self.dr)]))
        return True

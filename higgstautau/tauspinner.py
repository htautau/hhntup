"""
See instructions here:
https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/TauSpinner/trunk/README.rst

Based on original implementation by Daniele:
https://twiki.cern.ch/twiki/bin/view/Main/DanieleZanziTauSpinner
"""
from rootpy.tree.filtering import EventFilter
from TauSpinner import TauSpinnerTool, D3PD_MC

from . import log; log = log[__name__]


class EmbeddingTauSpinner(EventFilter):
    """
    This class applies Tau Spinner
    """
    def __init__(self, tree, passthrough=False, **kwargs):

        super(EmbeddingTauSpinner, self).__init__(
            passthrough=passthrough,
            **kwargs)

        if not passthrough:
            self.tree = tree
            self.spin_tool = TauSpinnerTool()
            self.spin_tool.initialize()

    def passes(self, event):
        # If passthrough is True (not embedding )then this
        # method is never called
        d3pd_mc = D3PD_MC()
        d3pd_mc.set_addresses(
           event['mc_n'],
           event.mc_pt,
           event.mc_eta,
           event.mc_phi,
           event.mc_m,
           event.mc_pdgId,
           event.mc_child_index)
        self.tree.embedding_spin_weight = self.spin_tool.get_spin_weight(d3pd_mc)
        return True

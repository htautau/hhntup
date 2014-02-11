"""
See instructions here:
https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/TauSpinner/trunk/README.rst

Based on original implementation by Daniele:
https://twiki.cern.ch/twiki/bin/view/Main/DanieleZanziTauSpinner
"""
from rootpy.tree.filtering import EventFilter
from . import log; log = log[__name__]

import ROOT


class EmbeddingTauSpinner(EventFilter):
    """
    This class applies Tau Spinner
    """
    def __init__(self, year, tree, passthrough=False, **kwargs):

        super(EmbeddingTauSpinner, self).__init__(
            passthrough=passthrough,
            **kwargs)

        if not passthrough:
            import TauSpinnerTool
            self.tree = tree
            self.d3pd_mc = ROOT.TauSpinnerHelpers.D3PD_MC()
            if year == 2012:
                self.tool = ROOT.TauSpinnerHelpers.HSG4.build8TeV()
            elif year == 2011:
                self.tool = ROOT.TauSpinnerHelpers.HSG4.build7TeV()
            else:
                raise ValueError(
                    "No TauSpinner config for year {0}".format(year))

    def passes(self, event):
        # If passthrough is True (not embedding) then this
        # method is never called
        self.d3pd_mc.set_addresses(
           event.mc_pt,
           event.mc_eta,
           event.mc_phi,
           event.mc_m,
           event.mc_pdgId,
           event.mc_status,
           event.mc_child_index)
        status = self.tool.read_event(self.d3pd_mc)
        self.tree.embedding_spin_weight = self.tool.get_spin_weight()
        return True

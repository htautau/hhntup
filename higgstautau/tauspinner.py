"""
See instructions here:
https://svnweb.cern.ch/trac/atlasphys/browser/Physics/Higgs/HSG4/software/common/TauSpinner/trunk/README.rst

https://twiki.cern.ch/twiki/bin/view/Main/DanieleZanziTauSpinner
"""
from atlastools import datasets
from rootpy.tree.filtering import EventFilter
from rootpy import ROOTError

# Load the library
from ROOT import gSystem
import os
gSystem.Load("./TauSpinner/lib/libTauSpinnerTool.so")
from ROOT import TauSpinnerTool, D3PD_MC
from . import log; log = log[__name__]

class EmbeddingTauSpinner(EventFilter):
    """
    This class applies Tau Spinner
    """
    def __init__(self, tree, datatype, year,
            verbose=False,
            passthrough=False,
            **kwargs):

        super(EmbeddingTauSpinner, self).__init__(
                passthrough=passthrough,
                **kwargs)

        if not passthrough:
            self.tree = tree
            self.year = year
            self.datatype = datatype
            self.isdata = datatype in (datasets.DATA, datasets.EMBED)
            self.verbose = verbose
            self.spin_tool = TauSpinnerTool()
            self.spin_tool.initialize()

    def passes(self, event):
        if not self.datatype is datasets.EMBED:
            # Not embedding -- do nothing!
            return True
        d3pd_mc = D3PD_MC()
        d3pd_mc.set_addresses(
           event['mc_n'],
           event.mc_pt,
           event.mc_eta,
           event.mc_phi,
           event.mc_m,
           event.mc_pdgId,
           event.mc_child_index)
        self.tree.spin_weight = self.spin_tool.get_spin_weight(d3pd_mc)
        return True

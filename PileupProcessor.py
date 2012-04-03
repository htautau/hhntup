"""
See here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ExtendedPileupReweighting#AnchorGen
"""

import ROOT
from atlastools.batch import ATLASStudent
from rootpy.tree import TreeChain

from higgstautau.pileup import PileupReweighting

ROOT.gErrorIgnoreLevel = ROOT.kFatal

# TODO: put this in HHSkim.py before next skimming

class PileupProcessor(ATLASStudent):

    def work(self):

        # initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                          files=self.files,
                          branches=[
                              'RunNumber',
                              'mc_channel_number',
                              'mc_event_weight',
                              'averageIntPerXing'
                              ])

        pileup_tool = PileupReweighting('pileup_weights')
        pileup_tool.UsePeriodConfig("MC11b")
        pileup_tool.Initialize()

        # entering the main event loop...
        for event in intree:

            pileup_tool.Fill(event.RunNumber,
                             event.mc_channel_number,
                             event.mc_event_weight,
                             event.averageIntPerXing)

        pileup_tool.WriteToFile()

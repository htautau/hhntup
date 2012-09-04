from rootpy.tree.filtering import EventFilter

from externaltools import PileupReweighting
from ROOT import Root
PileupReweighting = Root.TPileupReweighting


class PileupTemplates(EventFilter):

    def __init__(self, year, **kwargs):

        # initialize the pileup reweighting tool
        self.pileup_tool = PileupReweighting()
        if year == 2011:
            self.pileup_tool.UsePeriodConfig("MC11b")
        elif year == 2012:
            self.pileup_tool.UsePeriodConfig("MC12a")
        self.pileup_tool.Initialize()
        super(PileupTemplates, self).__init__(**kwargs)

    def passes(self, event):

        self.pileup_tool.Fill(
            event.RunNumber,
            event.mc_channel_number,
            event.mc_event_weight,
            event.averageIntPerXing)

    def finalize(self):

        # write the pileup reweighting file
        self.pileup_tool.WriteToFile()

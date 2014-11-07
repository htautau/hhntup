import ROOT
from rootpy.tree.filtering import EventFilter
from . import datasets
from .pileup import PILEUP_TOOLS
from .filters import BCH_TOOLS

RANDOMS = []

def get_random():
    # Note on use of ROOT random number generators:
    # TRandom and TRandom2 have many documented deficiencies.
    # TRandom3 is generally considered safely usable.
    # Also note that ROOT's gRandom calls TRandom3.
    random = ROOT.TRandom3()
    RANDOMS.append(random)
    return random


class RandomSeed(EventFilter):

    def __init__(self, datatype, **kwargs):
        super(RandomSeed, self).__init__(**kwargs)
        self.datatype = datatype

    def passes(self, event):
        if self.datatype in (datasets.DATA, datasets.EMBED):
            seed = int(event.EventInfo.runNumber() + event.EventInfo.eventNumber())
        else:
            seed = int(event.EventInfo.mcChannelNumber() + event.EventInfo.eventNumber())
        # METUtility uses gRandom
        ROOT.gRandom.SetSeed(seed)
        idx = 1
        for random in RANDOMS + BCH_TOOLS:
            random.SetSeed(seed + idx)
            idx += 1
        for pileup_tool in PILEUP_TOOLS:
            # same seeds for all pileup tools
            pileup_tool.SetRandomSeed(seed + idx)
        return True


class RandomRunNumber(EventFilter):

    def __init__(self, tree, datatype, pileup_tool, **kwargs):
        self.tree = tree
        self.pileup_tool = pileup_tool
        super(RandomRunNumber, self).__init__(**kwargs)
        if datatype in (datasets.MC, datasets.MCEMBED):
            self.passes = self.passes_mc
        else:
            self.passes = self.passes_data

    def passes_data(self, event):
        self.tree.RunNumber = event.EventInfo.runNumber()
        self.tree.lbn = event.EventInfo.lumiBlock()
        return True

    def passes_mc(self, event):
        # get random run number using the pileup tool
        # NEED TO UPDATE THE PILEUP TOOL FILE
        random_run = self.pileup_tool.GetRandomRunNumber(195847)
        # random_run = self.pileup_tool.GetRandomRunNumber(event.EventInfo.runNumber())
        self.tree.RunNumber = random_run
        self.tree.lbn = self.pileup_tool.GetRandomLumiBlockNumber(random_run)
        return True

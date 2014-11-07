# --
# November 6th, 2014: Not converted to XAOD yet
# --

from rootpy.tree.filtering import EventFilter

from externaltools import PileupReweighting
from ROOT import Root

from . import datasets
from . import log; log = log[__name__]
import os

# https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/InDetTrackingPerformanceGuidelines
PU_RESCALE = {
    2011: (0.97, 0.01),
    2012: (1.09, 0.04),
}

PILEUP_TOOLS = []


def get_pileup_reweighting_tool(year, use_defaults=True, systematic=None):
    # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ExtendedPileupReweighting
    # Initialize the pileup reweighting tool
    pileup_tool = Root.TPileupReweighting()
    if year == 2011:
        if use_defaults:
            pileup_tool.AddConfigFile(
                PileupReweighting.get_resource(
                    'mc11b_defaults.prw.root'))
        else:
            pileup_tool.AddConfigFile(
                'lumi/2011/'
                'TPileupReweighting.mc11.prw.root')
        lumicalc_file = 'lumi/2011/ilumicDDalc_histograms_None_178044-191933.root'
    elif year == 2012:
        if use_defaults:
            pileup_tool.AddConfigFile(
                PileupReweighting.get_resource(
                    'mc12ab_defaults.prw.root'))
        else:
            pileup_tool.AddConfigFile(
                'lumi/2012/'
                'TPileupReweighting.mc12.prw.root')
        lumicalc_file = 'lumi/2012/ilumicalc_histograms_None_200842-215643.root'
    else:
        raise ValueError(
            'No pileup reweighting defined for year %d' % year)
    rescale, rescale_error = PU_RESCALE[year]
    if systematic is None:
        pileup_tool.SetDataScaleFactors(1. / rescale)
    elif systematic == 'high':
        pileup_tool.SetDataScaleFactors(1. / (rescale + rescale_error))
    elif systematic == 'low':
        pileup_tool.SetDataScaleFactors(1. / (rescale - rescale_error))
    else:
        raise ValueError(
            "pileup systematic '{0}' not understood".format(systematic))
    pileup_tool.AddLumiCalcFile(lumicalc_file)
    # discard unrepresented data (with mu not simulated in MC)
    pileup_tool.SetUnrepresentedDataAction(2)
    pileup_tool.Initialize()
    # set the random seed used by the GetRandomRunNumber and
    # GetRandomPeriodNumber methods
    pileup_tool.SetRandomSeed(1777)
    # register
    PILEUP_TOOLS.append(pileup_tool)
    return pileup_tool


class PileupTemplates(EventFilter):

    def __init__(self, year, passthrough=False, **kwargs):
        if not passthrough:
            # initialize the pileup reweighting tool
            self.pileup_tool = Root.TPileupReweighting()
            if year == 2011:
                self.pileup_tool.UsePeriodConfig("MC11b")
            elif year == 2012:
                self.pileup_tool.UsePeriodConfig("MC12a")
            self.pileup_tool.Initialize()
        super(PileupTemplates, self).__init__(
            passthrough=passthrough,
            **kwargs)

    def passes(self, event):
        #pileup_chan107655_run195847
        self.pileup_tool.Fill(
            #     event.EventInfo.runNumber(),
            195847,
            event.EventInfo.mcChannelNumber(),
            event.EventInfo.mcEventWeight(),
            event.EventInfo.averageInteractionsPerCrossing())
        return True

    def finalize(self):
        if not self.passthrough:
            # write the pileup reweighting file
            self.pileup_tool.WriteToFile()


class PileupReweight(EventFilter):
    # XAOD MIGRATION: hard coding of the runnumber for now
    """
    Currently only implements hadhad reweighting
    """
    def __init__(self, year, tool, tool_high, tool_low,
                 tree, passthrough=False, **kwargs):
        if not passthrough:
            self.tree = tree
            self.tool = tool
            self.tool_high = tool_high
            self.tool_low = tool_low
        super(PileupReweight, self).__init__(
            passthrough=passthrough,
            **kwargs)

    def passes(self, event):
        # set the pileup weights
        self.tree.pileup_weight = self.tool.GetCombinedWeight(
            195847,
            # event.EventInfo.runNumber(),
            event.EventInfo.mcChannelNumber(),
            event.EventInfo.averageInteractionsPerCrossing())
        self.tree.pileup_weight_high = self.tool_high.GetCombinedWeight(
            195847,
            # event.EventInfo.runNumber(),
            event.EventInfo.mcChannelNumber(),
            event.EventInfo.averageInteractionsPerCrossing())
        self.tree.pileup_weight_low = self.tool_low.GetCombinedWeight(
            195847,
            # event.EventInfo.runNumber(),
            event.EventInfo.mcChannelNumber(),
            event.EventInfo.averageInteractionsPerCrossing())
        return True


class PileupScale(EventFilter):

    def __init__(self, tree, year, datatype, **kwargs):
        self.tree = tree
        self.scale = PU_RESCALE[year][0]
        super(PileupScale, self).__init__(**kwargs)
        if datatype in (datasets.DATA, datasets.EMBED):
            self.passes = self.passes_data
        elif datatype in (datasets.MC, datasets.MCEMBED):
            self.passes = self.passes_mc
        else:
            raise ValueError("no pileup scale defined for datatype %d" %
                datatype)

    def passes_data(self, event):
        self.tree.averageIntPerXing = event.EventInfo.averageInteractionsPerCrossing()
        self.tree.actualIntPerXing = event.EventInfo.actualInteractionsPerCrossing()
        return True

    def passes_mc(self, event):
        self.tree.averageIntPerXing = event.EventInfo.averageInteractionsPerCrossing() * self.scale
        self.tree.actualIntPerXing = event.EventInfo.actualInteractionsPerCrossing() * self.scale
        return True


class averageIntPerXingPatch(EventFilter):
    # NOT CONVERTED AND NOT NEEDED IN THE XAOD
    """
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/ExtendedPileupReweighting:

    NOTE (23/01/2013): A bug has been found in the d3pd making code, causing
    all MC12 samples to have a few of the averageIntPerXing values incorrectly set
    (some should be 0 but are set to 1). The bug does not affect data. To resolve
    this, when reading this branch, for both prw file generating and for when
    retrieving pileup weights, you should amend the value with the following line
    of code:

    averageIntPerXing = (isSimulation && lbn==1 && int(averageIntPerXing+0.5)==1) ? 0. : averageIntPerXing;
    """
    def passes(self, event):
        if event.lbn == 1 and int(event.averageIntPerXing + 0.5) == 1:
            event.averageIntPerXing = 0
        return True

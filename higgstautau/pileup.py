from rootpy.tree.filtering import EventFilter

from externaltools import PileupReweighting
from ROOT import Root


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

        self.pileup_tool.Fill(
            event.RunNumber,
            event.mc_channel_number,
            event.mc_event_weight,
            event.averageIntPerXing)
        return True

    def finalize(self):

        if not self.passthrough:
            # write the pileup reweighting file
            self.pileup_tool.WriteToFile()


class PileupReweight(EventFilter):

    def __init__(self, year, tree, passthrough=False, **kwargs):

        if not passthrough:
            self.tree = tree

            # Initialize the pileup reweighting tool
            self.pileup_tool = Root.TPileupReweighting()
            if year == 2011:
                self.pileup_tool.AddConfigFile(
                        PileupReweighting.get_resource(
                            'mc11b_defaults.prw.root'))
                self.pileup_tool.AddLumiCalcFile(
                        'lumi/2011/hadhad/'
                        'ilumicalc_histograms_None_178044-191933.root')
            elif year == 2012:
                self.pileup_tool.AddConfigFile(
                        PileupReweighting.get_resource(
                            'mc12a_defaults.prw.root'))
                self.pileup_tool.SetDataScaleFactors(1./1.11)
                self.pileup_tool.AddLumiCalcFile(
                        'lumi/2012/hadhad/'
                        'ilumicalc_histograms_None_200841-213250.root')
            else:
                raise ValueError('No pileup reweighting defined for year %d' %
                        year)
            # discard unrepresented data (with mu not simulated in MC)
            self.pileup_tool.SetUnrepresentedDataAction(2)
            self.pileup_tool.Initialize()

        super(PileupReweight, self).__init__(
                passthrough=passthrough,
                **kwargs)

    def passes(self, event):

        # set the event weight
        self.tree.pileup_weight = self.pileup_tool.GetCombinedWeight(
                event.RunNumber,
                event.mc_channel_number,
                event.averageIntPerXing)
        return True

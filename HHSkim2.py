import ROOT
import math

from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from atlastools.batch import ATLASStudent

from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel
from rootpy.types import *
from rootpy.io import open as ropen
from rootpy.plotting import Hist

from higgstautau.mixins import TauFourMomentum
from higgstautau.hadhad.filters import Triggers
import goodruns

#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class TriggerMatching(TreeModel):

    tau_trigger_match_index = ROOT.vector('int')
    tau_trigger_match_thresh = ROOT.vector('int')


class HHSkim2(ATLASStudent):

    def work(self):

        # merge TrigConfTrees
        metadirname = '%sMeta' % self.metadata.treename
        trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
        map(trigconfchain.Add, self.files)
        metadir = self.output.mkdir(metadirname)
        metadir.cd()
        trigconfchain.Merge(self.output, -1, 'fast keep')
        self.output.cd()

        # merge the cutflow hists from the first skim
        cutflow = None
        for file in self.files:
            with ropen(file) as f:
                if cutflow is None:
                    cutflow = file.cutflow.Clone()
                else:
                    cutflow += file.cutflow
        self.output.cd()
        cutflow.Write()

        if self.metadata.datatype == datasets.DATA:
            # merge GRL XML strings
            grls = []
            merged_grl = goodruns.GRL()
            for fname in self.files:
                merged_grl |= goodruns.GRL('%s:/Lumi/%s' % (fname, self.metadata.treename))
            lumi_dir = self.output.mkdir('Lumi')
            lumi_dir.cd()
            xml_string= ROOT.TObjString(merged_grl.str())
            xml_string.Write(self.metadata.treename)
            self.output.cd()

            # merge outtree_extras from the first skim
            extra_trees = ROOT.TChain(self.metadata.treename +
                                      '_failed_skim_after_trigger')
            map(extra_trees.Add, self.files)
            extra_trees.Merge(self.output, -1, 'fast keep')
            self.output.cd()

        onfilechange = []

        # initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                           files=self.files,
                           events=self.events,
                           onfilechange=onfilechange)

        Model = None
        if self.metadata.datatype == datasets.DATA:
            Model += TriggerMatching

        outtree = Tree(name=self.metadata.treename,
                       file=self.output,
                       model=Model)

        outtree.set_buffer(intree.buffer,
                           create_branches=True)

        # entering the main event loop...
        for event in intree:
            outtree.Fill()

        self.output.cd()

        # flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()

"""
Add instructions to run the skim here
"""

import ROOT
import math

from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from atlastools.batch import ATLASStudent

from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel, TreeBuffer
from rootpy.types import *
from rootpy.io import open as ropen
from rootpy.plotting import Hist
from rootpy.registry import lookup_demotion

from higgstautau.mixins import TauFourMomentum
from higgstautau.mixins import FourMomentum
from higgstautau.lephad.filters import AllmuTriggers, AlleTriggers, AllMCTriggers, \
                                       tau_skimselection, muon_skimselection, electron_skimselection
import goodruns


ROOT.gErrorIgnoreLevel = ROOT.kFatal


class SkimExtraModel(TreeModel):

    number_of_good_taus = IntCol()
    number_of_good_muons = IntCol()
    number_of_good_electrons = IntCol()


class LHSkim(ATLASStudent):

    def work(self):

        # Initialize the TreeChain of all input files
        intree = TreeChain(self.metadata.treename,
                          files=self.files,
                          events=self.events)

        outtree = Tree(name=self.metadata.treename,
                       file=self.output,
                       model=SkimExtraModel)
        outtree.set_buffer(intree.buffer, create_branches=True, visible=False)


        # Slimming
        blocks_to_remove= ['cl_*',
                           'trk_*',
                           'tau_otherTrk_*',
                           'jet_AntiKt6LCTopo_*',
                           'jet_AntiKt4TopoEM_*',
                           'mu_muid*',
                           'mu_calo*'
                           ]

        variables_to_keep = outtree.glob('*', prune=blocks_to_remove)
        outtree.activate(variables_to_keep, exclusive=True)


        # Merge TrigConfTrees
        metadirname = '%sMeta' % self.metadata.treename
        trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
        map(trigconfchain.Add, self.files)
        metadir = self.output.mkdir(metadirname)
        metadir.cd()
        trigconfchain.Merge(self.output, -1, 'fast keep')
        self.output.cd()


        # merge GRL XML strings
        if self.metadata.datatype == datasets.DATA:
            grls = []
            merged_grl = goodruns.GRL()
            for fname in self.files:
                merged_grl |= goodruns.GRL('%s:/Lumi/%s' % (fname, self.metadata.treename))
            lumi_dir = self.output.mkdir('Lumi')
            lumi_dir.cd()
            xml_string= ROOT.TObjString(merged_grl.str())
            xml_string.Write(self.metadata.treename)
            self.output.cd()


        # set the event filters
        trigger_filter = None
        if self.metadata.datatype == datasets.DATA:
            if self.metadata.title == datasets.MUON:
                trigger_filter = AllmuTriggers()
            if self.metadata.title == datasets.ELEC:
                trigger_filter = AlleTriggers()
        else:
            trigger_filter = AllMCTriggers()


        # Define collections for preselection
        intree.define_collection(name='taus', prefix='tau_', size='tau_n', mix=TauFourMomentum)
        intree.define_collection(name='muons', prefix='mu_staco_', size='mu_staco_n')
        intree.define_collection(name='electrons', prefix='el_', size='el_n')


        # Cut Flow counters
        nevents = 0
        nevents_mc_weight = 0
        nevents_passing_trigger = 0
        nevents_with_good_muons = 0
        nevents_with_good_electrons = 0
        nevents_with_good_taus = 0
        nevents_with_good_lephad = 0


        # Entering the main event loop...
        for event in intree:
            nevents += 1
            if self.metadata.datatype == datasets.MC:
                nevents_mc_weight += event.mc_event_weight
                
            if trigger_filter(event):
                nevents_passing_trigger +=1

                # Tau Skim Selection
                event.taus.select(lambda tau : tau_skimselection(tau))
                number_of_good_taus = len(event.taus)
                if number_of_good_taus > 0:
                    nevents_with_good_taus +=1

                # Muon Skim Selection
                event.muons.select(lambda mu : muon_skimselection(mu))
                number_of_good_muons = len(event.muons)
                if number_of_good_muons > 0:
                    nevents_with_good_muons +=1

                # Electron Skim Selection
                event.electrons.select(lambda e : electron_skimselection(e))
                number_of_good_electrons = len(event.electrons)
                if number_of_good_electrons > 0:
                    nevents_with_good_electrons +=1

                # Steer object selection for data and MC
                hasTau = number_of_good_taus > 0
                hasLeptonData = (number_of_good_muons > 0 and self.metadata.title == datasets.MUON) or \
                                (number_of_good_electrons > 0 and self.metadata.title == datasets.ELEC)
                hasLeptonMC = (number_of_good_muons > 0 or number_of_good_electrons > 0)
                isData = self.metadata.datatype == datasets.DATA
                isMC = self.metadata.datatype == datasets.MC

                # Fill the event
                if (isData and hasTau and hasLeptonData) or (isMC and hasTau and hasLeptonMC):
                    nevents_with_good_lephad +=1
                    outtree.number_of_good_taus = number_of_good_taus
                    outtree.number_of_good_muons = number_of_good_muons
                    outtree.number_of_good_electrons = number_of_good_electrons
                    outtree.Fill()


        self.output.cd()

        # Store the original number of events
        cutflow = Hist(7, 0, 7, name='cutflow', type='D')
        cutflow[0] = nevents
        cutflow[1] = nevents_mc_weights
        cutflow[2] = nevents_passing_trigger
        cutflow[3] = nevents_with_good_muons
        cutflow[4] = nevents_with_good_electrons
        cutflow[5] = nevents_with_good_taus
        cutflow[6] = nevents_with_good_lephad
        cutflow.Write()

        # Flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()

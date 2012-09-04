import ROOT
import math

from atlastools import utils
from atlastools import datasets
from atlastools.units import GeV
from atlastools.batch import ATLASStudent
from atlastools.filtering import GRLFilter

from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel
from rootpy.types import *
from rootpy.io import open as ropen

from higgstautau.mixins import TauFourMomentum, MCParticle
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.jetcalibration import JetCalibration
from higgstautau.patches import ElectronIDpatch, TauIDpatch
import higgstautau.skimming.hadhad as hhskimming
from hhskimming import branches as hhbranches
from hhskimming.models import *
from higgtautau.pileup import PileupTemplates, PileupReweight

import goodruns


#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class HHSkim(ATLASStudent):

    def work(self):

        if self.metadata.datatype != datasets.EMBED:
            # merge TrigConfTrees
            metadirname = '%sMeta' % self.metadata.treename
            trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
            map(trigconfchain.Add, self.files)
            metadir = self.output.mkdir(metadirname)
            metadir.cd()
            trigconfchain.Merge(self.output, -1, 'fast keep')
            self.output.cd()

        if self.metadata.datatype == datasets.DATA:
            # merge GRL XML strings
            grls = []
            merged_grl = goodruns.GRL()
            for fname in self.files:
                merged_grl |= goodruns.GRL(
                        '%s:/Lumi/%s' % (fname, self.metadata.treename))
            lumi_dir = self.output.mkdir('Lumi')
            lumi_dir.cd()
            xml_string= ROOT.TObjString(merged_grl.str())
            xml_string.Write(self.metadata.treename)
            self.output.cd()

        Model = SkimExtraModel
        if self.metadata.datatype == datasets.MC:
            Model += TriggerEmulation

        outtree = Tree(
                name=self.metadata.treename,
                file=self.output,
                model=Model)

        onfilechange = []
        count_funcs = {}

        if self.metadata.datatype == datasets.MC:
            # TODO get trigger emul tool here
            onfilechange.append(
                        (update_trigger_trees, (self, trigger_tool_wrapper,)))

            def mc_weight_count(event):
                return event.mc_event_weight

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        event_filters = EventFilterList([
            GRLFilter(
                self.grl,
                passthrough=(self.metadata.datatype != datasets.DATA
                             or self.metadata.year == 2012)),
            EmbeddingPileupPatch(
                passthrough=self.metadata.datatype != datasets.EMBED,
                count_funcs=count_funcs),
            PileupTemplates(
                passthrough=self.metadata.datatype != datasets.MC,
                count_funcs=count_funcs),
            #ExtraInfoTree(
            #   count_funcs=count_funcs)
            TauTriggerEmulation(
                year=self.metadata.year,
                outtree=outtree,
                passthrough=self.metadata.datatype != datasets.MC,
                count_funcs=count_funcs),
            Triggers(
                datatype=self.metadata.datatype,
                year=self.metadata.year,
                skim=True,
                count_funcs=count_funcs),
            PriVertex(
                count_funcs=count_funcs),
            LArError(
                count_funcs=count_funcs),
            # the BDT bits are broken in the p1130 production, correct them
            # DON'T FORGET TO REMOVE THIS WHEN SWITCHING TO A NEWER
            # PRODUCTION TAG!!!
            TauIDpatch(
                year=self.metadata.year,
                passthrough=self.metadata.year != 2012,
                count_funcs=count_funcs),
            # patch electron ID for 2012
            ElectronIDpatch(
                passthrough=self.metadata.year != 2012,
                count_funcs=count_funcs),
            # no need to recalibrate jets in 2012 (yet...)
            JetCalibration(
                datatype=self.metadata.datatype,
                year=self.metadata.year,
                verbose=VERBOSE,
                passthrough=self.metadata.year == 2012,
                count_funcs=count_funcs),
            LArHole(
                datatype=self.metadata.datatype,
                count_funcs=count_funcs),
            JetCleaning(
                datatype=self.metadata.datatype,
                year=self.metadata.year,
                count_funcs=count_funcs),
            ElectronVeto(
                count_funcs=count_funcs),
            MuonVeto(
                year=self.metadata.year,
                count_funcs=count_funcs),
            TauElectronVeto(2,
                count_funcs=count_funcs),
            TauMuonVeto(2,
                count_funcs=count_funcs),
            TauAuthor(2,
                count_funcs=count_funcs),
            TauHasTrack(2,
                count_funcs=count_funcs),
            TauPT(2,
                thresh=20 * GeV,
                count_funcs=count_funcs),
            TauEta(2,
                count_funcs=count_funcs),
            TauCrack(2,
                count_funcs=count_funcs),
            TauLArHole(2,
                count_funcs=count_funcs),
            TauID_BDTLoose_LLHLoose(2
                count_funcs=count_funcs),
            TauTriggerMatch(
                config=trigger_config,
                year=self.metadata.year,
                datatype=self.metadata.datatype,
                skim=True,
                tree=tree,
                passthrough=self.metadata.datatype == datasets.EMBED,
                count_funcs=count_funcs),
            PileupReweight(
                passthrough=self.metadata.datatype != datasets.MC,
                count_funcs=count_funcs),
        ])

        self.filters['event'] = event_filters

        # initialize the TreeChain of all input files
        chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                events=self.events,
                onfilechange=onfilechange,
                filters=event_filters)

        outtree.set_buffer(
                chain.buffer,
                ignore_branches=chain.glob(
                    hhbranches.REMOVE,
                    prune=hhbranches.KEEP),
                create_branches=True,
                visible=False)

        """
        if self.metadata.datatype == datasets.DATA:
            # outtree_extra holds info for events not included in the skim
            outtree_extra = Tree(
                    name=self.metadata.treename + '_failed_skim_after_trigger',
                    file=self.output,
                    model=SkimExtraModel + SkimExtraTauPtModel)

            extra_variables = [
                'trig_EF_tau_pt',
                'actualIntPerXing',
                'averageIntPerXing',
                'MET_RefFinal_BDTMedium_phi',
                'MET_RefFinal_BDTMedium_et'
                'MET_RefFinal_BDTMedium_sumet',
                'EventNumber',
                'RunNumber',
                'lbn'
            ]

            if self.metadata.year % 1000 == 11:
                extra_variables += Triggers.triggers_11
            elif self.metadata.year % 1000 == 12:
                extra_variables += Triggers.triggers_12
            else:
                raise ValueError("No triggers defined for year %d" % year)

            outtree_extra.set_buffer(
                    chain.buffer,
                    branches=extra_variables,
                    create_branches=True,
                    visible=False)
        else:
            # write out some branches for all events
            # that failed the skim before trigger only for MC
            # Used to get the total pileup reweighting sum
            if self.metadata.datatype == datasets.MC:
                outtree_extra = Tree(
                        name=self.metadata.treename + '_failed_skim_before_trigger',
                        file=self.output)
            else: #embedding
                outtree_extra = Tree(
                        name=self.metadata.treename + '_failed_skim_before_selection',
                        file=self.output)
            extra_variables = [
                'actualIntPerXing',
                'averageIntPerXing',
                'RunNumber',
            ]
            outtree_extra.set_buffer(
                    chain.buffer,
                    branches=extra_variables,
                    create_branches=True, visible=False)
        """

        # define tree collections
        chain.define_collection(
                name="taus",
                prefix="tau_",
                size="tau_n",
                mix=TauFourMomentum)
        chain.define_collection(
                name="taus_EF",
                prefix="trig_EF_tau_",
                size="trig_EF_tau_n",
                mix=TauFourMomentum)
        # jet_* etc. is AntiKt4LCTopo_* in tau-perf D3PDs
        chain.define_collection(
                name="jets",
                prefix="jet_",
                size="jet_n",
                mix=FourMomentum)
        chain.define_collection(
                name="jets_EM",
                prefix="jet_AntiKt4TopoEM_",
                size="jet_AntiKt4TopoEM_n",
                mix=FourMomentum)
        chain.define_collection(
                name="truetaus",
                prefix="trueTau_",
                size="trueTau_n",
                mix=MCTauFourMomentum)
        chain.define_collection(
                name="mc",
                prefix="mc_",
                size="mc_n",
                mix=MCParticle)
        chain.define_collection(
                name="muons",
                prefix="mu_staco_",
                size="mu_staco_n")
        chain.define_collection(
                name="electrons",
                prefix="el_",
                size="el_n")
        chain.define_collection(
                name="vertices",
                prefix="vxp_",
                size="vxp_n")

        # entering the main event loop...
        for event in chain:

            event.vertices.select(vertex_selection)
            outtree.number_of_good_vertices = len(event.vertices)
            outtree.Fill()
            """
            elif self.metadata.datatype == datasets.DATA:
                outtree_extra.number_of_good_vertices = number_of_good_vertices
                if event.taus:
                    # There can be at most one good tau if this event failed the skim
                    outtree_extra.tau_pt = event.taus[0].pt
                else:
                    outtree_extra.tau_pt = -1111.
                outtree_extra.Fill()
            elif self.metadata.datatype == datasets.EMBED:
                outtree_extra.Fill()
            """

        self.output.cd()

        """
        if self.metadata.datatype == datasets.MC:
            # store the original weighted number of events
            cutflow = Hist(2, 0, 2, name='cutflow', type='D')
            cutflow[1] = nevents_mc_weight
        else:
            cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        # store the original number of events
        cutflow[0] = nevents
        cutflow.Write()
        """

        # flush any baskets remaining in memory to disk
        outtree.FlushBaskets()
        outtree.Write()
        """
        outtree_extra.FlushBaskets()
        outtree_extra.Write()
        """

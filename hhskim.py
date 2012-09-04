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

from higgstautau.mixins import *
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau.embedding import EmbeddingPileupPatch
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.emulation import TauTriggerEmulation
from higgstautau.jetcalibration import JetCalibration
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.skimming.hadhad import branches as hhbranches
from higgstautau.skimming.hadhad.models import *
from higgstautau.pileup import PileupTemplates, PileupReweight

import goodruns


#ROOT.gErrorIgnoreLevel = ROOT.kFatal
VALIDATE = False
VERBOSE = False


class hhskim(ATLASStudent):

    def work(self):

        datatype = self.metadata.datatype
        year = self.metadata.year

        if datatype != datasets.EMBED:
            # merge TrigConfTrees
            metadirname = '%sMeta' % self.metadata.treename
            trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
            map(trigconfchain.Add, self.files)
            metadir = self.output.mkdir(metadirname)
            metadir.cd()
            trigconfchain.Merge(self.output, -1, 'fast keep')
            self.output.cd()

        if datatype == datasets.DATA:
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

        Model = SkimModel + TriggerMatching
        if datatype == datasets.MC:
            if year == 2011:
                Model += TriggerEmulation11

        tree = Tree(
                name=self.metadata.treename,
                file=self.output,
                model=Model)

        onfilechange = []
        count_funcs = {}

        if datatype == datasets.MC:

            def mc_weight_count(event):
                return event.mc_event_weight

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        trigger_emulation = TauTriggerEmulation(
                year=year,
                tree=tree,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),

        if datatype == datasets.MC:
            onfilechange.append(
                (trigger_emulation.update_trigger_trees,
                    (self, trigger_emulation)))

        trigger_config = None

        if datatype != datasets.EMBED:
            # trigger config tool to read trigger info in the ntuples
            trigger_config = get_trigger_config()

            # update the trigger config maps on every file change
            onfilechange.append((update_trigger_config, (trigger_config,)))

        event_filters = EventFilterList([
            GRLFilter(
                self.grl,
                passthrough=(datatype != datasets.DATA
                             or year == 2012),
                count_funcs=count_funcs),
            EmbeddingPileupPatch(
                passthrough=datatype != datasets.EMBED,
                count_funcs=count_funcs),
            PileupTemplates(
                year=year,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            #ExtraInfoTree(
            #   count_funcs=count_funcs)
            Triggers(
                datatype=datatype,
                year=year,
                count_funcs=count_funcs),
            PriVertex(
                count_funcs=count_funcs),
            LArError(
                count_funcs=count_funcs),
            # the BDT bits are broken in the p1130 production, correct them
            # DON'T FORGET TO REMOVE THIS WHEN SWITCHING TO A NEWER
            # PRODUCTION TAG!!!
            TauIDpatch(
                year=year,
                passthrough=year != 2012,
                count_funcs=count_funcs),
            # patch electron ID for 2012
            ElectronIDpatch(
                passthrough=year != 2012,
                count_funcs=count_funcs),
            # no need to recalibrate jets in 2012 (yet...)
            JetCalibration(
                datatype=datatype,
                year=year,
                verbose=VERBOSE,
                passthrough=year == 2012,
                count_funcs=count_funcs),
            LArHole(
                datatype=datatype,
                count_funcs=count_funcs),
            JetCleaning(
                datatype=datatype,
                year=year,
                count_funcs=count_funcs),
            ElectronVeto(
                count_funcs=count_funcs),
            MuonVeto(
                year=year,
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
            TauID_BDTLoose_LLHLoose(2,
                count_funcs=count_funcs),
            TauTriggerMatch(
                config=trigger_config,
                year=year,
                datatype=datatype,
                skim=True,
                tree=tree,
                passthrough=datatype == datasets.EMBED,
                count_funcs=count_funcs),
            PileupReweight(
                year=year,
                tree=tree,
                passthrough=datatype != datasets.MC,
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

        tree.set_buffer(
                chain.buffer,
                ignore_branches=chain.glob(
                    hhbranches.REMOVE,
                    prune=hhbranches.KEEP),
                create_branches=True,
                visible=False)

        """
        if datatype == datasets.DATA:
            # tree_extra holds info for events not included in the skim
            tree_extra = Tree(
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

            if year % 1000 == 11:
                extra_variables += Triggers.triggers_11
            elif year % 1000 == 12:
                extra_variables += Triggers.triggers_12
            else:
                raise ValueError("No triggers defined for year %d" % year)

            tree_extra.set_buffer(
                    chain.buffer,
                    branches=extra_variables,
                    create_branches=True,
                    visible=False)
        else:
            # write out some branches for all events
            # that failed the skim before trigger only for MC
            # Used to get the total pileup reweighting sum
            if datatype == datasets.MC:
                tree_extra = Tree(
                        name=self.metadata.treename + '_failed_skim_before_trigger',
                        file=self.output)
            else: #embedding
                tree_extra = Tree(
                        name=self.metadata.treename + '_failed_skim_before_selection',
                        file=self.output)
            extra_variables = [
                'actualIntPerXing',
                'averageIntPerXing',
                'RunNumber',
            ]
            tree_extra.set_buffer(
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

        if VALIDATE: # only validate on a single data run or MC channel
            chain.GetEntry(0)
            if datatype == datasets.MC:
                validate_log = open('skim2_validate_mc_%d.txt' %
                        chain.mc_channel_number, 'w')
            elif datatype == datasets.DATA:
                validate_log = open('skim2_validate_data_%d.txt' %
                        chain.RunNumber, 'w')
            else:
                validate_log = open('skim2_validate_embedded_%d.txt' %
                        chain.RunNumber, 'w')

        # entering the main event loop...
        for event in chain:

            event.vertices.select(vertex_selection)
            tree.number_of_good_vertices = len(event.vertices)

            if datatype == datasets.EMBED:
                # select two leading taus by pT
                event.taus.sort(key=lambda tau: tau.pt, reverse=True)
                event.taus.slice(stop=2)

            assert len(event.taus) == 2

            selected_idx = [tau.index for tau in event.taus]
            selected_idx.sort()

            if VALIDATE:
                if datatype == datasets.MC:
                    print >> validate_log, event.mc_channel_number,
                print >> validate_log, event.RunNumber, event.EventNumber,
                print >> validate_log, "%.4f" % tree.pileup_weight,
                for idx in selected_idx:
                    print >> validate_log, idx, tree.tau_trigger_match_thresh[idx],
                print >> validate_log

            tree.tau_selected.clear()
            for i in xrange(event.tau_n):
                if i in selected_idx:
                    tree.tau_selected.push_back(True)
                else:
                    tree.tau_selected.push_back(False)

            tree.Fill()

            """
            elif datatype == datasets.DATA:
                tree_extra.number_of_good_vertices = number_of_good_vertices
                if event.taus:
                    # There can be at most one good tau if this event failed the skim
                    tree_extra.tau_pt = event.taus[0].pt
                else:
                    tree_extra.tau_pt = -1111.
                tree_extra.Fill()
            elif datatype == datasets.EMBED:
                tree_extra.Fill()
            """

        self.output.cd()

        """
        if datatype == datasets.MC:
            # store the original weighted number of events
            cutflow = Hist(2, 0, 2, name='cutflow', type='D')
            cutflow[1] = nevents_mc_weight
        else:
            cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        # store the original number of events
        cutflow[0] = nevents
        cutflow.Write()
        """

        if VALIDATE:
            validate_log.close()
            # sort the output by event number for MC
            # sort +2 -3 -n skim2_validate_mc_125205.txt -o skim2_validate_mc_125205.txt

        # flush any baskets remaining in memory to disk
        tree.FlushBaskets()
        tree.Write()
        """
        tree_extra.FlushBaskets()
        tree_extra.Write()
        """

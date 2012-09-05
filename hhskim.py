import ROOT

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
from higgstautau.trigger.emulation import TauTriggerEmulation, update_trigger_trees
from higgstautau.jetcalibration import JetCalibration
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.skimming.hadhad import branches as hhbranches
from higgstautau.skimming.hadhad.models import *
from higgstautau.pileup import PileupTemplates, PileupReweight
from higgstautau.hadhad.collections import define_collections

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

        tree = Tree(
                name=self.metadata.treename,
                file=self.output,
                model=SkimModel + TriggerMatching)

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
                count_funcs=count_funcs)

        if datatype == datasets.MC:
            onfilechange.append(
                (update_trigger_trees, (self, trigger_emulation,)))

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
            trigger_emulation,
            Triggers(
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
                filters=event_filters,
                verbose=True)

        define_collections(chain)

        tree.set_buffer(
                chain.buffer,
                ignore_branches=chain.glob(
                    hhbranches.REMOVE,
                    prune=hhbranches.KEEP),
                create_branches=True,
                transfer_objects=True)

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

        self.output.cd()

        if VALIDATE:
            validate_log.close()
            # sort the output by event number for MC
            # sort +2 -3 -n skim2_validate_mc_125205.txt -o skim2_validate_mc_125205.txt

        # flush any baskets remaining in memory to disk
        tree.FlushBaskets()
        tree.Write()

import ROOT

import os
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
from rootpy.extern.argparse import ArgumentParser

from higgstautau.mixins import *
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau import mass
from higgstautau.embedding import EmbeddingPileupPatch, EmbeddingIsolation
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.emulation import TauTriggerEmulation, update_trigger_trees
from higgstautau.trigger.matching import TauTriggerMatchIndex, TauTriggerMatchThreshold
from higgstautau.trigger.efficiency import TauTriggerEfficiency
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.skimming.hadhad import branches as hhbranches
from higgstautau.skimming.hadhad.models import *
from higgstautau.hadhad.models import EmbeddingBlock
from higgstautau.pileup import PileupTemplates, PileupReweight
from higgstautau.hadhad.objects import define_objects
from higgstautau.corrections import reweight_ggf

import goodruns


#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class hhskim(ATLASStudent):

    def __init__(self, options, **kwargs):

        super(hhskim, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-terms', default=None)
        parser.add_argument('--no-trigger', action='store_true', default=False)
        parser.add_argument('--no-grl', action='store_true', default=False)
        parser.add_argument('--student-verbose', action='store_true', default=False)
        parser.add_argument('--validate', action='store_true', default=False)
        self.args = parser.parse_args(options)
        if self.args.syst_terms is not None:
            self.args.syst_terms = set([
                eval('Systematics.%s' % term) for term in
                self.args.syst_terms.split(',')])

    def work(self):

        datatype = self.metadata.datatype
        year = self.metadata.year
        no_trigger = self.args.no_trigger
        no_grl = self.args.no_grl
        verbose = self.args.student_verbose
        validate = self.args.validate

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

        # define the model of the output tree
        Model = SkimModel + TriggerMatching + MassModel + TauCorrections

        if datatype == datasets.EMBED:
            # add embedding systematics branches
            Model += EmbeddingBlock

        self.output.cd()

        # create the output tree
        tree = Tree(
                name=self.metadata.treename,
                model=Model)

        onfilechange = []
        count_funcs = {}

        if datatype in (datasets.MC, datasets.EMBED):

            def mc_weight_count(event):
                return event.mc_event_weight

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        trigger_emulation = TauTriggerEmulation(
                year=year,
                passthrough=no_trigger or datatype != datasets.MC or year > 2011,
                count_funcs=count_funcs)

        if not trigger_emulation.passthrough:
            onfilechange.append(
                (update_trigger_trees, (self, trigger_emulation,)))

        trigger_config = None

        if datatype != datasets.EMBED:
            # trigger config tool to read trigger info in the ntuples
            trigger_config = get_trigger_config()

            # update the trigger config maps on every file change
            onfilechange.append((update_trigger_config, (trigger_config,)))

        # define the list of event filters
        event_filters = EventFilterList([
            CoreFlags(
                count_funcs=count_funcs), #TODO move this below GRL
            GRLFilter(
                self.grl,
                passthrough=(no_grl or datatype != datasets.DATA),
                count_funcs=count_funcs),
            EmbeddingPileupPatch(
                passthrough=year > 2011 or datatype != datasets.EMBED,
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
                passthrough=no_trigger or datatype == datasets.EMBED,
                count_funcs=count_funcs),
            PriVertex(
                count_funcs=count_funcs),
            LArError(
                count_funcs=count_funcs),
            TileError(
                count_funcs=count_funcs),
            JetCalibration(
                datatype=datatype,
                year=year,
                verbose=verbose,
                count_funcs=count_funcs),
            # PUT THE SYSTEMATICS "FILTER" BEFORE
            # ANY FILTERS THAT REFER TO OBJECTS
            # BUT AFTER CALIBRATIONS
            Systematics(
                terms=self.args.syst_terms,
                year=year,
                datatype=datatype,
                verbose=verbose,
                count_funcs=count_funcs),
            # the BDT bits are broken in the p1130 production, correct them
            # DON'T FORGET TO REMOVE THIS WHEN SWITCHING TO A NEWER
            # PRODUCTION TAG!!!
            #TauIDpatch(
            #    year=year,
            #    count_funcs=count_funcs),
            # patch electron ID for 2012
            #ElectronIDpatch(
            #    passthrough=year != 2012,
            #    count_funcs=count_funcs),
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
            TauTriggerMatchIndex(
                config=trigger_config,
                year=year,
                datatype=datatype,
                passthrough=no_trigger or datatype == datasets.EMBED,
                count_funcs=count_funcs),
            # select two leading taus at this point
            # 25/35 for data
            # 20/30 for MC for TES
            TauLeadSublead(
                lead=35 * GeV if datatype == datasets.DATA else 30 * GeV,
                sublead=25 * GeV if datatype == datasets.DATA else 20 * GeV,
                count_funcs=count_funcs),
            TaudR(3.2,
                count_funcs=count_funcs),
            TauTriggerMatchThreshold(
                passthrough=no_trigger or datatype == datasets.EMBED,
                count_funcs=count_funcs),
            TauTriggerEfficiency(
                year=year,
                datatype=datatype,
                tes_systematic=self.args.syst_terms and (
                    Systematics.TES_TERMS & self.args.syst_terms),
                passthrough=no_trigger or datatype == datasets.DATA,
                count_funcs=count_funcs),
            PileupReweight(
                year=year,
                tree=tree,
                use_defaults=True,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            TruthMatching(
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            EfficiencyScaleFactors(
                year=year,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            FakeRateScaleFactors(
                year=year,
                passthrough=no_trigger or datatype != datasets.MC
                            or year == 2012, # wait for new tool
                count_funcs=count_funcs),
            ggFReweighting(
                dsname=os.getenv('INPUT_DATASET_NAME', ''),
                tree=tree,
                # no ggf reweighting for 2012 MC
                passthrough=datatype != datasets.MC or year != 2011,
                count_funcs=count_funcs),
            TauTrackRecounting(
                year=year,
                count_funcs=count_funcs),
            EmbeddingIsolation(
                tree=tree,
                passthrough=year < 2012 or datatype != datasets.EMBED,
                count_funcs=count_funcs),
        ])

        # set the event filters
        self.filters['event'] = event_filters

        # initialize the TreeChain of all input files
        chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                events=self.events,
                onfilechange=onfilechange,
                filters=event_filters,
                cache=True,
                cache_size=50000000,
                learn_entries=100)

        define_objects(chain, year, skim=True)

        # include the branches in the input chain in the output tree
        # set branches to be removed in ignore_branches
        tree.set_buffer(
                chain._buffer,
                ignore_branches=chain.glob(
                    hhbranches.REMOVE,
                    exclude=hhbranches.KEEP),
                create_branches=True,
                transfer_objects=True,
                visible=False)

        tree.define_object(name='tau', prefix='tau_')

        if validate: # only validate on a single data run or MC channel
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

        # create MMC object
        mmc = mass.MMC(year=year, channel='hh')

        self.output.cd()
        # entering the main event loop...
        for event in chain:

            assert len(event.taus) == 2

            tree.number_of_good_vertices = len(event.vertices)

            tau1, tau2 = event.taus
            selected_idx = [tau.index for tau in event.taus]
            selected_idx.sort()

            if validate:
                if datatype == datasets.MC:
                    print >> validate_log, event.mc_channel_number,
                print >> validate_log, event.RunNumber, event.EventNumber,
                print >> validate_log, "%.4f" % tree.pileup_weight,
                for idx in selected_idx:
                    print >> validate_log, idx, tree.tau_trigger_match_thresh[idx],
                print >> validate_log

            # missing mass
            METx = event.MET.etx
            METy = event.MET.ety
            MET = event.MET.et
            sumET = event.MET.sumet

            mmc_mass, mmc_resonance, mmc_met = mmc.mass(
                    tau1, tau2, METx, METy, sumET)

            tree.tau_MMC_mass = mmc_mass
            tree.tau_MMC_resonance.set_from(mmc_resonance)
            if mmc_mass > 0:
                tree.tau_MMC_resonance_pt = mmc_resonance.Pt()
            tree.tau_MMC_MET = mmc_met.Mod()
            tree.tau_MMC_MET_x = mmc_met.X()
            tree.tau_MMC_MET_y = mmc_met.Y()
            tree.tau_MMC_MET_phi = math.pi - mmc_met.Phi()

            # collinear mass
            collin_mass, tau1_x, tau2_x = mass.collinearmass(
                    tau1, tau2, METx, METy)
            tree.tau_collinear_mass = collin_mass
            tau1.collinear_momentum_fraction = tau1_x
            tau2.collinear_momentum_fraction = tau2_x

            # visible mass
            tree.tau_visible_mass = (tau1.fourvect + tau2.fourvect).M()

            if datatype == datasets.MC and year == 2011:
                tree.ggf_weight = reweight_ggf(event, self.metadata.name)

            tree.tau_selected.clear()
            tree.tau_collinear_momentum_fraction.clear()

            SkimModel.reset(tree)
            TriggerMatching.reset(tree)
            TauCorrections.reset(tree)

            # set the skim-defined variables in the output tree
            for i in xrange(event.tau_n):
                tau = event.taus.getitem(i)

                tree.tau_selected.push_back(i in selected_idx)
                tree.tau_collinear_momentum_fraction.push_back(
                    tau.collinear_momentum_fraction)

                SkimModel.set(tree, tau)
                TriggerMatching.set(tree, tau)
                TauCorrections.set(tree, tau)

            # fill the output tree
            tree.Fill()

        self.output.cd()

        if validate:
            validate_log.close()
            # sort the output by event number for MC
            # sort +2 -3 -n skim2_validate_mc_125205.txt -o skim2_validate_mc_125205.txt

        # flush any baskets remaining in memory to disk
        tree.FlushBaskets()
        tree.Write()

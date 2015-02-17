# Author: Noel Dawe
# Modified by: Quentin Buat (xAOD migration)

import ROOT

import os
import math
import subprocess

import goodruns
# import externaltools

# rootpy imports
from rootpy.plotting import Hist
from rootpy.tree.filtering import EventFilter, EventFilterList
from rootpy.tree import Tree, TreeChain, TreeModel, TreeBuffer
from rootpy.extern.argparse import ArgumentParser
from rootpy.io import root_open
from rootpy import stl, asrootpy
from rootpy.vector import LorentzVector

# local xaod imports
from xaod.xaodtree import xAODTree
# local higgstautau imports
from higgstautau import eventshapes, utils, datasets
from higgstautau.batch import ATLASStudent
from higgstautau.units import GeV
from higgstautau.filters import *
from higgstautau.hadhad.objects import define_objects
from higgstautau.hadhad.models import *
from higgstautau.hadhad.filters import *
from higgstautau import mass
from higgstautau.mass import is_MET_bisecting
# from higgstautau.embedding import *
# from higgstautau.systematics import Systematics
# from higgstautau.met import METRecalculation
from higgstautau.jetcalibration import JetCalibration, JetResolution
# from higgstautau.tauspinner import EmbeddingTauSpinner
# from higgstautau.trigger import update_trigger_config, get_trigger_config
# from higgstautau.trigger.efficiency import TauTriggerEfficiency
# from higgstautau.trigger.emulation import (
#     TauTriggerEmulation, update_trigger_trees)
from higgstautau.pileup import PileupScale, PileupReweight_xAOD
from higgstautau.rand import RandomRunNumber, RandomSeed
from higgstautau import log; log = log[__name__]


class hhskim(ATLASStudent):

    def __init__(self, options, **kwargs):
        super(hhskim, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--local', action='store_true', default=False)
        parser.add_argument('--syst-terms', default=None)
        parser.add_argument('--student-verbose', action='store_true', default=False)
        parser.add_argument('--student-very-verbose', action='store_true', default=False)
        parser.add_argument('--redo-selection', action='store_true', default=False)
        parser.add_argument('--nominal-values', action='store_true', default=False)
        args = parser.parse_args(options)
        self.args = args
        if args.syst_terms is not None:
            args.syst_terms = set([
                eval('Systematics.%s' % term) for term in
                args.syst_terms.split(',')])
        if args.local:
            def merge(inputs, output, metadata):
                # merge output trees
                root_output = output + '.root'
                log.info("merging output trees")
                subprocess.call(['hadd', root_output] + inputs)
                if metadata.datatype == datasets.DATA:
                    # merge GRLs
                    log.info("merging GRL fragments")
                    grl = goodruns.GRL()
                    for input in inputs:
                        grl |= goodruns.GRL('%s:/lumi' % input)
                    grl.save('%s:/lumi' % root_output)

            hhskim.merge = staticmethod(merge)

    def work(self):
        # get argument values
        local = self.args.local
        syst_terms = self.args.syst_terms
        datatype = self.metadata.datatype
        year = self.metadata.year
        verbose = self.args.student_verbose
        very_verbose = self.args.student_very_verbose
        redo_selection = self.args.redo_selection
        nominal_values = self.args.nominal_values

        # get the dataset name
        dsname = os.getenv('INPUT_DATASET_NAME', None)
        if dsname is None:
            # attempt to guess dsname from dirname
            if self.files:
                dsname = os.path.basename(os.path.dirname(self.files[0]))

        # is this a signal sample?
        # if so we will also keep some truth information in the output below
        is_signal = datatype == datasets.MC and (
            '_VBFH' in dsname or
            '_ggH' in dsname or
            '_ZH' in dsname or
            '_WH' in dsname or
            '_ttH' in dsname)
        log.info("DATASET: {0}".format(dsname))
        log.info("IS SIGNAL: {0}".format(is_signal))

        # is this an inclusive signal sample for overlap studies?
        is_inclusive_signal = is_signal and '_inclusive' in dsname

        # is this a BCH-fixed sample? (temporary)
        is_bch_sample = 'r5470_r4540_p1344' in dsname
        if is_bch_sample:
            log.warning("this is a BCH-fixed r5470 sample")

        # onfilechange will contain a list of functions to be called as the
        # chain rolls over to each new file
        onfilechange = []
        count_funcs = {}

        if datatype != datasets.DATA:
            # count the weighted number of events
            if local:
                def mc_weight_count(event):
                    return event.hh_mc_weight
            else:
                def mc_weight_count(event):
                    return event.TruthEvent[0].weights()[0]

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        if local:
            # local means running on the skims, the output of this script
            # running on the grid
            if datatype == datasets.DATA:
                # merge the GRL fragments
                merged_grl = goodruns.GRL()

                def update_grl(student, grl, name, file, tree):
                    grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

                onfilechange.append((update_grl, (self, merged_grl,)))

            if datatype == datasets.DATA:
                merged_cutflow = Hist(1, 0, 1, name='cutflow', type='D')
            else:
                merged_cutflow = Hist(2, 0, 2, name='cutflow', type='D')

            def update_cutflow(student, cutflow, name, file, tree):
                # record a cut-flow
                year = student.metadata.year
                datatype = student.metadata.datatype
                cutflow[1].value += file.cutflow_event[1].value
                if datatype != datasets.DATA:
                    cutflow[2].value += file.cutflow_event_mc_weight[1].value

            onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        else:

            # NEED TO BE CONVERTED TO XAOD
            # if datatype not in (datasets.EMBED, datasets.MCEMBED):
            #     # merge TrigConfTrees
            #     metadirname = '%sMeta' % self.metadata.treename
            #     trigconfchain = ROOT.TChain('%s/TrigConfTree' % metadirname)
            #     map(trigconfchain.Add, self.files)
            #     metadir = self.output.mkdir(metadirname)
            #     metadir.cd()
            #     trigconfchain.Merge(self.output, -1, 'fast keep')
            #     self.output.cd()

            if datatype == datasets.DATA:
                # merge GRL XML strings
                merged_grl = goodruns.GRL()
            #     for fname in self.files:
            #         with root_open(fname) as f:
            #             for key in f.Lumi.keys():
            #                 merged_grl |= goodruns.GRL(
            #                     str(key.ReadObj().GetString()),
            #                     from_string=True)
            #     lumi_dir = self.output.mkdir('Lumi')
            #     lumi_dir.cd()
            #     xml_string= ROOT.TObjString(merged_grl.str())
            #     xml_string.Write(self.metadata.treename)
            #     self.output.cd()

        self.output.cd()

        # create the output tree
        model = get_model(datatype, dsname,
                          prefix=None if local else 'hh_',
                          is_inclusive_signal=is_inclusive_signal)
        log.info("Output Model:\n\n{0}\n\n".format(model))
        outtree = Tree(name=self.metadata.treename,
                       model=model)

        if local:
            tree = outtree
        else:
            tree = outtree.define_object(name='tree', prefix='hh_')

        #tree.define_object(name='tau', prefix='tau_')
        tree.define_object(name='tau1', prefix='tau1_')
        tree.define_object(name='tau2', prefix='tau2_')
        tree.define_object(name='truetau1', prefix='truetau1_')
        tree.define_object(name='truetau2', prefix='truetau2_')
        tree.define_object(name='jet1', prefix='jet1_')
        tree.define_object(name='jet2', prefix='jet2_')
        tree.define_object(name='jet3', prefix='jet3_')

        mmc_objects = [
            tree.define_object(name='mmc0', prefix='mmc0_'),
            tree.define_object(name='mmc1', prefix='mmc1_'),
            tree.define_object(name='mmc2', prefix='mmc2_'),
        ]

        for mmc_obj in mmc_objects:
            mmc_obj.define_object(name='resonance', prefix='resonance_')

        # NEED TO BE CONVERTED TO XAOD
        # trigger_emulation = TauTriggerEmulation(
        #     year=year,
        #     passthrough=local or datatype != datasets.MC or year > 2011,
        #     count_funcs=count_funcs)

        # if not trigger_emulation.passthrough:
        #     onfilechange.append(
        #         (update_trigger_trees, (self, trigger_emulation,)))

        # trigger_config = None

        # if datatype not in (datasets.EMBED, datasets.MCEMBED):
        #     # trigger config tool to read trigger info in the ntuples
        #     trigger_config = get_trigger_config()
        #     # update the trigger config maps on every file change
        #     onfilechange.append((update_trigger_config, (trigger_config,)))

        # define the list of event filters
        if local and syst_terms is None and not redo_selection:
            event_filters = None
        else:
            tau_ntrack_recounted_use_ntup = False
            if year > 2011:
                # peek at first tree to determine if the extended number of
                # tracks is already stored
                with root_open(self.files[0]) as test_file:
                    test_tree = test_file.Get(self.metadata.treename)
                    tau_ntrack_recounted_use_ntup = (
                        'tau_out_track_n_extended' in test_tree)

            log.info(self.grl)
            event_filters = EventFilterList([
                GRLFilter(
                    self.grl,
                    passthrough=(
                        local or (
                            datatype not in (datasets.DATA, datasets.EMBED))),
                    count_funcs=count_funcs),
                CoreFlags(
                    passthrough=local,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # EmbeddingPileupPatch(
                #     passthrough=(
                #         local or year > 2011 or datatype != datasets.EMBED),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD (not a priority)
                # PileupTemplates(
                #     year=year,
                #     passthrough=(
                #         local or is_bch_sample or datatype not in (
                #             datasets.MC, datasets.MCEMBED)),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # RandomSeed(
                #     datatype=datatype,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # BCHSampleRunNumber(
                #     passthrough=not is_bch_sample,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # RandomRunNumber(
                #     tree=tree,
                #     datatype=datatype,
                #     pileup_tool=pileup_tool,
                #     passthrough=local,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # trigger_emulation,
                # NEED TO BE CONVERTED TO XAOD
                # Triggers(
                #     year=year,
                #     tree=tree,
                #     datatype=datatype,
                #     passthrough=datatype in (datasets.EMBED, datasets.MCEMBED),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                PileupReweight_xAOD(
                        tree=tree,
                        passthrough=(local or (
                            datatype not in (datasets.MC, datasets.MCEMBED))),
                        count_funcs=count_funcs),
                PriVertex(
                    passthrough=local,
                    count_funcs=count_funcs),
                LArError(
                    passthrough=local,
                    count_funcs=count_funcs),
                TileError(
                    passthrough=local,
                    count_funcs=count_funcs),
                TileTrips(
                    passthrough=(
                        local or datatype in (datasets.MC, datasets.MCEMBED)),
                    count_funcs=count_funcs),
                JetCalibration(
                        datatype=datatype,
                        passthrough=local,
                        count_funcs=count_funcs),
                JetResolution(
                        passthrough=(local or (
                                datatype not in (datasets.MC, datasets.MCEMBED))),
                        count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # # in situ TES shift for 2012 data
                # TauEnergyShift(
                #     passthrough=(
                #         local or datatype != datasets.DATA
                #         or year < 2012 or nominal_values),
                #     count_funcs=count_funcs),
                # # truth matching must come before systematics due to
                # # TES_TRUE/FAKE
                # NEED TO BE CONVERTED TO XAOD
                TrueTauSelection(
                        passthrough=datatype == datasets.DATA,
                        count_funcs=count_funcs),
                TruthMatching(
                    passthrough=datatype == datasets.DATA,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                NvtxJets(
                    tree=tree,
                    count_funcs=count_funcs),
                # # PUT THE SYSTEMATICS "FILTER" BEFORE
                # # ANY FILTERS THAT REFER TO OBJECTS
                # # BUT AFTER CALIBRATIONS
                # # Systematics must also come before anything that refers to
                # # thing.fourvect since fourvect is cached!
                # NEED TO BE CONVERTED TO XAOD
                # Systematics(
                #     terms=syst_terms,
                #     year=year,
                #     datatype=datatype,
                #     tree=tree,
                #     verbose=verbose,
                #     passthrough=not syst_terms,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # JetIsPileup(
                #     passthrough=(
                #         local or year < 2012 or
                #         datatype not in (datasets.MC, datasets.MCEMBED)),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # LArHole(
                #     tree=tree,
                #     passthrough=year > 2011,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                JetCleaning(
                    datatype=datatype,
                    year=year,
                    count_funcs=count_funcs),
                # Need to check the electron ID and OQ
                ElectronVeto(
                        el_sel='Medium',
                        count_funcs=count_funcs),
                MuonVeto(
                    count_funcs=count_funcs),
                TauPT(2,
                    thresh=20 * GeV,
                    count_funcs=count_funcs),
                TauHasTrack(2,
                    count_funcs=count_funcs),
                TauEta(2,
                    count_funcs=count_funcs),
                TauElectronVeto(2,
                    count_funcs=count_funcs),
                TauMuonVeto(2,
                    count_funcs=count_funcs),
                TauCrack(2,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # TauLArHole(2,
                #     tree=tree,
                #     passthrough=year > 2011,
                #     count_funcs=count_funcs),
                # # before selecting the leading and subleading taus
                # # be sure to only consider good candidates
                TauIDMedium(2,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # but not used by default
                # #TauTriggerMatchIndex(
                # #    config=trigger_config,
                # #    year=year,
                # #    datatype=datatype,
                # #    passthrough=datatype == datasets.EMBED,
                # #    count_funcs=count_funcs),
                # Select two leading taus at this point
                # 25 and 35 for data
                # 20 and 30 for MC to leave room for TES uncertainty
                TauLeadSublead(
                    lead=(
                        35 * GeV if datatype == datasets.DATA or local
                        else 30 * GeV),
                    sublead=(
                        25 * GeV if datatype == datasets.DATA or local
                        else 20 * GeV),
                    count_funcs=count_funcs),
                # taus are sorted (in decreasing order) by pT from here on
                TauIDSelection(
                    count_funcs=count_funcs),
                TaudR(3.2,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # but not used by default
                # #TauTriggerMatchThreshold(
                # #    datatype=datatype,
                # #    tree=tree,
                # #    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # TauTriggerEfficiency(
                #     year=year,
                #     datatype=datatype,
                #     tree=tree,
                #     tes_systematic=self.args.syst_terms and (
                #         Systematics.TES_TERMS & self.args.syst_terms),
                #     passthrough=datatype == datasets.DATA,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                PileupScale(
                    tree=tree,
                    year=year,
                    datatype=datatype,
                    passthrough=local,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                TauIDScaleFactors(
                    year=year,
                    passthrough=datatype == datasets.DATA,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # TauFakeRateScaleFactors(
                #     year=year,
                #     datatype=datatype,
                #     tree=tree,
                #     tes_up=(self.args.syst_terms is not None and
                #         (Systematics.TES_FAKE_TOTAL_UP in self.args.syst_terms or
                #          Systematics.TES_FAKE_FINAL_UP in self.args.syst_terms)),
                #     tes_down=(self.args.syst_terms is not None and
                #         (Systematics.TES_FAKE_TOTAL_DOWN in self.args.syst_terms or
                #          Systematics.TES_FAKE_FINAL_DOWN in self.args.syst_terms)),
                #     passthrough=datatype in (datasets.DATA, datasets.EMBED),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                HiggsPT(
                    year=year,
                    tree=tree,
                    passthrough=not is_signal or local,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # TauTrackRecounting(
                #     year=year,
                #     use_ntup_value=tau_ntrack_recounted_use_ntup,
                #     passthrough=local,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # MCWeight(
                #     datatype=datatype,
                #     tree=tree,
                #     passthrough=local or datatype == datasets.DATA,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # EmbeddingIsolation(
                #     tree=tree,
                #     passthrough=(
                #         local or year < 2012 or
                #         datatype not in (datasets.EMBED, datasets.MCEMBED)),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # EmbeddingCorrections(
                #     tree=tree,
                #     year=year,
                #     passthrough=(
                #         local or
                #         datatype not in (datasets.EMBED, datasets.MCEMBED)),
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # EmbeddingTauSpinner(
                #     year=year,
                #     tree=tree,
                #     passthrough=(
                #         local or datatype not in (
                #             datasets.EMBED, datasets.MCEMBED)),
                #     count_funcs=count_funcs),
                # # put MET recalculation after tau selection but before tau-jet
                # # overlap removal and jet selection because of the RefAntiTau
                # # MET correction
                # NEED TO BE CONVERTED TO XAOD
                # METRecalculation(
                #     terms=syst_terms,
                #     year=year,
                #     tree=tree,
                #     refantitau=not nominal_values,
                #     verbose=verbose,
                #     very_verbose=very_verbose,
                #     count_funcs=count_funcs),
                TauJetOverlapRemoval(
                    count_funcs=count_funcs),
                JetPreselection(
                    count_funcs=count_funcs),
                NonIsolatedJet(
                    tree=tree,
                    count_funcs=count_funcs),
                JetSelection(
                    year=year,
                    count_funcs=count_funcs),
                RecoJetTrueTauMatching(
                    passthrough=datatype == datasets.DATA or local,
                    count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # BCHCleaning(
                #     tree=tree,
                #     passthrough=year == 2011 or local,
                #     datatype=datatype,
                #     count_funcs=count_funcs),
                # NEED TO BE CONVERTED TO XAOD
                # ClassifyInclusiveHiggsSample(
                #     tree=tree,
                #     passthrough=not is_inclusive_signal,
                #     count_funcs=count_funcs),
            ])

            # set the event filters
            self.filters['event'] = event_filters

        hh_buffer = TreeBuffer()
        if local:
            chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                # ignore_branches=ignore_branches,
                events=self.events,
                onfilechange=onfilechange,
                filters=event_filters,
                cache=True,
                cache_size=50000000,
                learn_entries=100)
            buffer = TreeBuffer()
            for name, value in chain._buffer.items():
                if name.startswith('hh_'):
                    hh_buffer[name[3:]] = value
                elif name in copied:
                    buffer[name] = value
            outtree.set_buffer(
                hh_buffer,
                create_branches=False,
                visible=True)
            outtree.set_buffer(
                buffer,
                create_branches=True,
                visible=False)


        else:

            root_chain = ROOT.TChain(self.metadata.treename)
            for f in self.files:
                log.info(f)
                root_chain.Add(f)
            
            # if len(self.files) != 1:
            #     raise RuntimeError('lenght of files has to be 1 for now (no xAOD chaining available)')
            # self.files = self.files[0]
            # root_chain = ROOT.TFile(self.files)

            chain = xAODTree(root_chain, filters=event_filters, events=self.events)
            define_objects(chain, datatype=datatype)
            outtree.set_buffer(
                hh_buffer,
                create_branches=True,
                visible=False)

            # create the MMC
            mmc = mass.MMC(year=year)

        # report which packages have been loaded
        # externaltools.report()

        self.output.cd()

        # The main event loop
        # the event filters above are automatically run for each event and only
        # the surviving events are looped on
        for event in chain:

            if local and syst_terms is None and not redo_selection:
                outtree.Fill()
                continue
            
            # sort taus and jets in decreasing order by pT
            event.taus.sort(key=lambda tau: tau.pt(), reverse=True)
            event.jets.sort(key=lambda jet: jet.pt(), reverse=True)

            # tau1 is the leading tau
            # tau2 is the subleading tau
            tau1, tau2 = event.taus
            tau1.fourvect = asrootpy(tau1.p4())
            tau2.fourvect = asrootpy(tau2.p4())

            beta_taus = (tau1.fourvect + tau2.fourvect).BoostVector()
            tau1.fourvect_boosted = LorentzVector()
            tau1.fourvect_boosted.copy_from(tau1.fourvect)
            tau1.fourvect_boosted.Boost(beta_taus * -1)
            
            tau2.fourvect_boosted = LorentzVector()
            tau2.fourvect_boosted.copy_from(tau2.fourvect)
            tau2.fourvect_boosted.Boost(beta_taus * -1)

            jets = list(event.jets)
            for jet in jets:
                jet.fourvect = asrootpy(jet.p4())

            jet1, jet2, jet3 = None, None, None
            beta = None
            if len(jets) >= 2:
                jet1, jet2 = jets[:2]

                # determine boost of system
                # determine jet CoM frame
                beta = (jet1.fourvect + jet2.fourvect).BoostVector()
                tree.jet_beta.copy_from(beta)

                jet1.fourvect_boosted = LorentzVector()
                jet1.fourvect_boosted.copy_from(jet1.fourvect)
                jet1.fourvect_boosted.Boost(beta * -1)

                jet2.fourvect_boosted = LorentzVector()
                jet2.fourvect_boosted.copy_from(jet2.fourvect)
                jet2.fourvect_boosted.Boost(beta * -1)

                tau1.min_dr_jet = min(
                    tau1.fourvect.DeltaR(jet1.fourvect),
                    tau1.fourvect.DeltaR(jet2.fourvect))
                tau2.min_dr_jet = min(
                    tau2.fourvect.DeltaR(jet1.fourvect),
                    tau2.fourvect.DeltaR(jet2.fourvect))

                # tau centrality (degree to which they are between the two jets)
                tau1.centrality = eventshapes.eta_centrality(
                    tau1.fourvect.Eta(),
                    jet1.fourvect.Eta(),
                    jet2.fourvect.Eta())

                tau2.centrality = eventshapes.eta_centrality(
                    tau2.fourvect.Eta(),
                    jet1.fourvect.Eta(),
                    jet2.fourvect.Eta())

                # boosted tau centrality
                tau1.centrality_boosted = eventshapes.eta_centrality(
                    tau1.fourvect_boosted.Eta(),
                    jet1.fourvect_boosted.Eta(),
                    jet2.fourvect_boosted.Eta())

                tau2.centrality_boosted = eventshapes.eta_centrality(
                    tau2.fourvect_boosted.Eta(),
                    jet1.fourvect_boosted.Eta(),
                    jet2.fourvect_boosted.Eta())

                # 3rd leading jet
                if len(jets) >= 3:
                    jet3 = jets[2]
                    jet3.fourvect_boosted = LorentzVector()
                    jet3.fourvect_boosted.copy_from(jet3.fourvect)
                    jet3.fourvect_boosted.Boost(beta * -1)

            elif len(jets) == 1:
                jet1 = jets[0]

                tau1.min_dr_jet = tau1.fourvect.DeltaR(jet1.fourvect)
                tau2.min_dr_jet = tau2.fourvect.DeltaR(jet1.fourvect)

            RecoJetBlock.set(tree, jet1, jet2, jet3, local=local)

            # mass of ditau + leading jet system
            if jet1 is not None:
                tree.mass_tau1_tau2_jet1 = (
                    tau1.fourvect + tau2.fourvect + jet1.fourvect).M()

            #####################################
            # number of tracks from PV minus taus
            #####################################
            ntrack_pv = 0
            ntrack_nontau_pv = 0
            for vxp in event.vertices:
                # primary vertex
                if vxp.vertexType() == 1:
                    ntrack_pv = vxp.nTrackParticles()
                    ntrack_nontau_pv = ntrack_pv - tau1.nTracks() - tau2.nTracks()
                    break
            tree.ntrack_pv = ntrack_pv
            tree.ntrack_nontau_pv = ntrack_nontau_pv

            #########################
            # MET variables
            #########################
            MET = event.MET[0]
            METx = MET.mpx()
            METy = MET.mpy()
            METet = MET.met()
            MET_vect = Vector2(METx, METy)
            MET_4vect = LorentzVector()
            MET_4vect.SetPxPyPzE(METx, METy, 0., METet)
            MET_4vect_boosted = LorentzVector()
            MET_4vect_boosted.copy_from(MET_4vect)
            if beta is not None:
                MET_4vect_boosted.Boost(beta * -1)

            tree.MET_et = METet
            tree.MET_etx = METx
            tree.MET_ety = METy
            tree.MET_phi = MET.phi()
            dPhi_tau1_tau2 = abs(tau1.fourvect.DeltaPhi(tau2.fourvect))
            dPhi_tau1_MET = abs(tau1.fourvect.DeltaPhi(MET_4vect))
            dPhi_tau2_MET = abs(tau2.fourvect.DeltaPhi(MET_4vect))
            tree.dPhi_tau1_tau2 = dPhi_tau1_tau2
            tree.dPhi_tau1_MET = dPhi_tau1_MET
            tree.dPhi_tau2_MET = dPhi_tau2_MET
            tree.dPhi_min_tau_MET = min(dPhi_tau1_MET, dPhi_tau2_MET)
            tree.MET_bisecting = is_MET_bisecting(
                dPhi_tau1_tau2,
                dPhi_tau1_MET,
                dPhi_tau2_MET)

            sumET = MET.sumet()
            tree.MET_sumet = sumET
            if sumET != 0:
                tree.MET_sig = ((2. * METet / GeV) /
                    (utils.sign(sumET) * sqrt(abs(sumET / GeV))))
            else:
                tree.MET_sig = -1.

            tree.MET_centrality = eventshapes.phi_centrality(
                tau1.fourvect,
                tau2.fourvect,
                MET_vect)
            tree.MET_centrality_boosted = eventshapes.phi_centrality(
                tau1.fourvect_boosted,
                tau2.fourvect_boosted,
                MET_4vect_boosted)

            tree.number_of_good_vertices = len(event.vertices)

            ##########################
            # Jet and sum pt variables
            ##########################
            tree.numJets = len(event.jets)

            # sum pT with only the two leading jets
            tree.sum_pt = sum(
                [tau1.pt(), tau2.pt()] +
                [jet.pt() for jet in jets[:2]])

            # sum pT with all selected jets
            tree.sum_pt_full = sum(
                [tau1.pt(), tau2.pt()] +
                [jet.pt() for jet in jets])

            # vector sum pT with two leading jets and MET
            tree.vector_sum_pt = sum(
                [tau1.fourvect, tau2.fourvect] +
                [jet.fourvect for jet in jets[:2]] +
                [MET_4vect]).Pt()

            # vector sum pT with all selected jets and MET
            tree.vector_sum_pt_full = sum(
                [tau1.fourvect, tau2.fourvect] +
                [jet.fourvect for jet in jets] +
                [MET_4vect]).Pt()

            # resonance pT
            tree.resonance_pt = sum(
                [tau1.fourvect, tau2.fourvect, MET_4vect]).Pt()

            # #############################
            # # tau <-> vertex association
            # #############################
            tree.tau_same_vertex = (
                tau1.vertex() == tau2.vertex())

            tau1.vertex_prob = ROOT.TMath.Prob(
                tau1.vertex().chiSquared(),
                int(tau1.vertex().numberDoF()))

            tau2.vertex_prob = ROOT.TMath.Prob(
                tau2.vertex().chiSquared(),
                int(tau2.vertex().numberDoF()))

            # ##########################
            # # MMC Mass
            # ##########################
            mmc_result = mmc.mass(
                tau1, tau2,
                METx, METy, sumET,
                njets=len(event.jets))

            for mmc_method, mmc_object in enumerate(mmc_objects):
                mmc_mass, mmc_resonance, mmc_met = mmc_result[mmc_method]
                if verbose:
                    log.info("MMC (method %d): %f" % (mmc_method, mmc_mass))

                mmc_object.mass = mmc_mass
                mmc_object.MET_et = mmc_met.Mod()
                mmc_object.MET_etx = mmc_met.X()
                mmc_object.MET_ety = mmc_met.Y()
                mmc_object.MET_phi = math.pi - mmc_met.Phi()
                if mmc_mass > 0:
                    FourMomentum.set(mmc_object.resonance, mmc_resonance)

            # ############################
            # # collinear and visible mass
            # ############################
            vis_mass, collin_mass, tau1_x, tau2_x = mass.collinearmass(
                tau1, tau2, METx, METy)

            tree.mass_vis_tau1_tau2 = vis_mass
            tree.mass_collinear_tau1_tau2 = collin_mass
            tau1.collinear_momentum_fraction = tau1_x
            tau2.collinear_momentum_fraction = tau2_x

            # # Fill the tau block
            # # This must come after the RecoJetBlock is filled since
            # # that sets the jet_beta for boosting the taus
            RecoTauBlock.set(event, tree, datatype, tau1, tau2, local=local)

            # NEED TO BE CONVERTED TO XAOD
            if datatype != datasets.DATA:
                TrueTauBlock.set(tree, tau1, tau2)
            # fill the output tree
            outtree.Fill(reset=True)

        # externaltools.report()

        # flush any baskets remaining in memory to disk
        self.output.cd()
        outtree.FlushBaskets()
        outtree.Write()

        if local:
            if datatype == datasets.DATA:
                xml_string = ROOT.TObjString(merged_grl.str())
                xml_string.Write('lumi')
            merged_cutflow.Write()

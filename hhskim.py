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
from rootpy.extern.argparse import ArgumentParser
from rootpy.io import root_open

from higgstautau import hepmc
from higgstautau import tautools
from higgstautau import eventshapes
from higgstautau.mixins import *
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau import mass
from higgstautau.mass import is_MET_bisecting
from higgstautau.overlap import TauJetOverlapRemoval
from higgstautau.embedding import EmbeddingPileupPatch, EmbeddingIsolation
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.emulation import (TauTriggerEmulation,
                                           update_trigger_trees)
from higgstautau.trigger.matching import (TauTriggerMatchIndex,
                                          TauTriggerMatchThreshold)
from higgstautau.trigger.efficiency import TauTriggerEfficiency
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.hadhad import branches as hhbranches
from higgstautau.hadhad.models import *
from higgstautau.pileup import (PileupTemplates, PileupReweight,
                                get_pileup_reweighting_tool,
                                averageIntPerXingPatch, PileupScale)
from higgstautau.hadhad.objects import define_objects
from higgstautau.corrections import reweight_ggf
from higgstautau import log; log = log[__name__]

import goodruns


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

        # get pileup reweighting tool
        pileup_tool = get_pileup_reweighting_tool(
            year=year,
            use_defaults=True)

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

        self.output.cd()

        # create the output tree
        tree = Tree(
            name=self.metadata.treename,
            model=get_model(datatype, self.metadata.name))

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
                passthrough=(
                    no_grl or datatype not in (datasets.DATA, datasets.EMBED)),
                count_funcs=count_funcs),
            EmbeddingPileupPatch(
                passthrough=year > 2011 or datatype != datasets.EMBED,
                count_funcs=count_funcs),
            averageIntPerXingPatch(
                passthrough=year < 2012 or datatype != datasets.MC,
                count_funcs=count_funcs),
            PileupTemplates(
                year=year,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            trigger_emulation,
            Triggers(
                year=year,
                passthrough=no_trigger or datatype == datasets.EMBED,
                count_funcs=count_funcs),
            PileupReweight(
                tool=pileup_tool,
                tree=tree,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            RandomRunNumber(
                tree=tree,
                datatype=datatype,
                pileup_tool=pileup_tool,
                count_funcs=count_funcs),
            PriVertex(
                count_funcs=count_funcs),
            LArError(
                count_funcs=count_funcs),
            TileError(
                count_funcs=count_funcs),
            TileTrips(
                passthrough=year < 2012,
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
            TauAuthor(2,
                count_funcs=count_funcs),
            TauCrack(2,
                count_funcs=count_funcs),
            TauLArHole(2,
                count_funcs=count_funcs),
            TauID_SkimLoose(2,
                year=year,
                count_funcs=count_funcs),
            TauTriggerMatchIndex(
                config=trigger_config,
                year=year,
                datatype=datatype,
                passthrough=no_trigger or datatype == datasets.EMBED,
                count_funcs=count_funcs),
            # Select two leading taus at this point
            # 25 and 35 for data
            # 20 and 30 for MC for TES uncertainty
            TauLeadSublead(
                lead=35 * GeV if datatype == datasets.DATA else 30 * GeV,
                sublead=25 * GeV if datatype == datasets.DATA else 20 * GeV,
                count_funcs=count_funcs),
            # apply this selection here since skim has lower threshold for data
            TauLeadSublead(
                lead=35 * GeV,
                sublead=25 * GeV,
                count_funcs=count_funcs),
            # taus are sorted (in decreasing order) by pT from here on
            TauIDSelection(
                year=year,
                tree=tree,
                count_funcs=count_funcs),
            TaudR(3.2,
                count_funcs=count_funcs),
            TruthMatching(
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            TauTriggerMatchThreshold(
                datatype=datatype,
                tree=tree,
                passthrough=no_trigger,
                count_funcs=count_funcs),
            TauTriggerEfficiency(
                year=year,
                datatype=datatype,
                tree=tree,
                tes_systematic=self.args.syst_terms and (
                    Systematics.TES_TERMS & self.args.syst_terms),
                passthrough=no_trigger or datatype == datasets.DATA,
                count_funcs=count_funcs),
            PileupScale(
                tree=tree,
                year=year,
                datatype=datatype,
                count_funcs=count_funcs),
            EfficiencyScaleFactors(
                year=year,
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            FakeRateScaleFactors(
                year=year,
                datatype=datatype,
                tes_up_systematic=(self.args.syst_terms and
                    (Systematics.TES_UP in self.args.syst_terms)),
                tes_down_systematic=(self.args.syst_terms and
                    (Systematics.TES_DOWN in self.args.syst_terms)),
                passthrough=no_trigger or datatype == datasets.DATA,
                count_funcs=count_funcs),
            ggFReweighting(
                dsname=os.getenv('INPUT_DATASET_NAME', ''),
                tree=tree,
                # no ggf reweighting for 2012 MC
                passthrough=datatype != datasets.MC or year != 2011,
                count_funcs=count_funcs),
            TauTrackRecounting(
                year=year,
                datatype=datatype,
                count_funcs=count_funcs),
            MCWeight(
                datatype=datatype,
                tree=tree,
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            EmbeddingIsolation(
                tree=tree,
                passthrough=year < 2012 or datatype != datasets.EMBED,
                count_funcs=count_funcs),
            EmbeddingCorrections(
                tree=tree,
                passthrough=year < 2012 or datatype != datasets.EMBED,
                count_funcs=count_funcs),
            TauJetOverlapRemoval(
                count_funcs=count_funcs),
            JetPreselection(
                passthrough=year < 2012,
                count_funcs=count_funcs),
            NonIsolatedJet(
                tree=tree,
                passthrough=year < 2012,
                count_funcs=count_funcs),
            JetSelection(
                year=year,
                count_funcs=count_funcs),
        ])

        # set the event filters
        self.filters['event'] = event_filters

        # peek at first tree to determine which branches to exclude
        with root_open(self.files[0]) as test_file:
            test_tree = test_file.Get(self.metadata.treename)
            ignore_branches = test_tree.glob(
                hhbranches.REMOVE,
                exclude=hhbranches.KEEP)

        # initialize the TreeChain of all input files
        chain = TreeChain(
            self.metadata.treename,
            files=self.files,
            ignore_branches=ignore_branches,
            events=self.events,
            onfilechange=onfilechange,
            filters=event_filters,
            cache=True,
            cache_size=50000000,
            learn_entries=100)

        # include the branches in the input chain in the output tree
        # set branches to be removed in ignore_branches
        tree.set_buffer(
            chain._buffer,
            ignore_branches=ignore_branches,
            create_branches=True,
            ignore_duplicates=True,
            transfer_objects=True,
            visible=False)

        if validate: # only validate on a single data run or MC channel
            chain.GetEntry(0)
            if datatype == datasets.MC:
                validate_log = open('hhskim_validate_mc_%d.txt' %
                    chain.mc_channel_number, 'w')
            elif datatype == datasets.DATA:
                validate_log = open('hhskim_validate_data_%d.txt' %
                    chain.RunNumber, 'w')
            else:
                validate_log = open('hhskim_validate_embedded_%d.txt' %
                    chain.RunNumber, 'w')

        # define tree objects
        define_objects(chain, year)

        tree.define_object(name='tau', prefix='tau_')
        tree.define_object(name='tau1', prefix='tau1_')
        tree.define_object(name='tau2', prefix='tau2_')
        tree.define_object(name='jet1', prefix='jet1_')
        tree.define_object(name='jet2', prefix='jet2_')

        mmc_objects = [
            tree.define_object(name='mmc0', prefix='mmc0_'),
            tree.define_object(name='mmc1', prefix='mmc1_'),
            tree.define_object(name='mmc2', prefix='mmc2_'),
        ]

        # create the MMC
        mmc = mass.MMC(year=year)

        self.output.cd()

        #####################
        # The main event loop
        #####################
        for event in chain:

            # sort taus and jets in decreasing order by pT
            event.taus.sort(key=lambda tau: tau.pt, reverse=True)
            event.jets.sort(key=lambda jet: jet.pt, reverse=True)

            # tau1 is the leading tau
            # tau2 is the subleading tau
            tau1, tau2 = event.taus
            jets = list(event.jets)

            if len(jets) >= 2:
                jet1, jet2 = jets[:2]

                # determine boost of system
                # determine jet CoM frame
                beta = (jet1.fourvect + jet2.fourvect).BoostVector()
                tree.jet_beta.set_from(beta)

                jet1.fourvect_boosted.set_from(jet1.fourvect)
                jet2.fourvect_boosted.set_from(jet2.fourvect)
                jet1.fourvect_boosted.Boost(beta * -1)
                jet2.fourvect_boosted.Boost(beta * -1)

                tau1.fourvect_boosted.set_from(tau1.fourvect)
                tau2.fourvect_boosted.set_from(tau2.fourvect)
                tau1.fourvect_boosted.Boost(beta * -1)
                tau2.fourvect_boosted.Boost(beta * -1)

                RecoJetBlock.set(tree, jet1, jet2)

                tau1.min_dr_jet = min(
                    tau1.fourvect.DeltaR(jet1.fourvect),
                    tau1.fourvect.DeltaR(jet2.fourvect))
                tau2.min_dr_jet = min(
                    tau2.fourvect.DeltaR(jet1.fourvect),
                    tau2.fourvect.DeltaR(jet2.fourvect))

                sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                    [tau1.fourvect,
                     tau2.fourvect,
                     jet1.fourvect,
                     jet2.fourvect])

                # sphericity
                tree.sphericity = sphericity
                # aplanarity
                tree.aplanarity = aplanarity

                sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                    [tau1.fourvect_boosted,
                     tau2.fourvect_boosted,
                     jet1.fourvect_boosted,
                     jet2.fourvect_boosted])

                # sphericity
                tree.sphericity_boosted = sphericity
                # aplanarity
                tree.aplanarity_boosted = aplanarity

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

            elif len(jets) >= 1:
                jet1 = jets[0]
                RecoJetBlock.set(tree, jet1)

                tau1.min_dr_jet = tau1.fourvect.DeltaR(jet1.fourvect)
                tau2.min_dr_jet = tau2.fourvect.DeltaR(jet1.fourvect)

                sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                    [tau1.fourvect,
                     tau2.fourvect,
                     jet1.fourvect])

                # sphericity
                tree.sphericity = sphericity
                # aplanarity
                tree.aplanarity = aplanarity

            #####################################
            # number of tracks from PV minus taus
            #####################################
            ntrack_pv = 0
            ntrack_nontau_pv = 0
            for vxp in event.vertices:
                # primary vertex
                if vxp.type == 1:
                    ntrack_pv = vxp.nTracks
                    ntrack_nontau_pv = ntrack_pv - tau1.numTrack - tau2.numTrack
                    break
            tree.ntrack_pv = ntrack_pv
            tree.ntrack_nontau_pv = ntrack_nontau_pv

            #########################
            # MET variables
            #########################
            METx = event.MET.etx
            METy = event.MET.ety
            MET = event.MET.et
            MET_vect = Vector2(METx, METy)
            MET_4vect = LorentzVector()
            MET_4vect.SetPxPyPzE(METx, METy, 0., MET)
            # TODO: save this in the output ntuple

            tree.MET = MET
            tree.MET_x = METx
            tree.MET_y = METy
            tree.MET_phi = event.MET.phi
            tree.MET_vec.set_from(MET_vect)
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

            sumET = event.MET.sumet
            tree.sumET = sumET
            if sumET != 0:
                tree.MET_sig = ((2. * MET / GeV) /
                        (utils.sign(sumET) * sqrt(abs(sumET / GeV))))
            else:
                tree.MET_sig = -1.

            tree.MET_centrality = eventshapes.phi_centrality(
                    tau1.fourvect,
                    tau2.fourvect,
                    MET_vect)

            tree.number_of_good_vertices = len(event.vertices)
            tau1, tau2 = event.taus

            selected_idx = [tau.index for tau in event.taus]
            selected_idx.sort()

            ##########################
            # Jet and sum pt variables
            ##########################
            tree.numJets = len(event.jets)

            # sum pT with only the two leading jets
            tree.sum_pt = sum(
                [tau1.pt, tau2.pt] +
                [jet.pt for jet in jets[:2]])

            # sum pT with all selected jets
            tree.sum_pt_full = sum(
                [tau1.pt, tau2.pt] +
                [jet.pt for jet in jets])

            # vector sum pT with two leading jets and MET
            tree.vector_sum_pt = sum(
                [tau1.fourvect, tau2.fourvect] +
                [jet.fourvect for jet in jets[:2]] +
                [MET_4vect]).Pt()

            # resonance pT
            tree.resonance_pt = sum(
                [tau1.fourvect, tau2.fourvect, MET_4vect]).Pt()

            #############################
            # tau <-> vertex association
            #############################
            tree.tau_same_vertex = (
                tau1.privtx_x == tau2.privtx_x and
                tau1.privtx_y == tau2.privtx_y and
                tau1.privtx_z == tau2.privtx_z)

            tau1.vertex_prob = ROOT.TMath.Prob(
                tau1.privtx_chiSquared,
                int(tau1.privtx_numberDoF))

            tau2.vertex_prob = ROOT.TMath.Prob(
                tau2.privtx_chiSquared,
                int(tau2.privtx_numberDoF))

            ##########################
            # MMC Mass
            ##########################
            METx = event.MET.etx
            METy = event.MET.ety
            MET = event.MET.et
            sumET = event.MET.sumet

            mmc_result = mmc.mass(
                tau1, tau2,
                METx, METy, sumET,
                njets=len(event.jets))

            for mmc_method, mmc_object in enumerate(mmc_objects):
                mmc_mass, mmc_resonance, mmc_met = mmc_result[mmc_method]
                if verbose:
                    log.info("MMC (method %d): %f" % (mmc_method, mmc_mass))

                mmc_object.mass = mmc_mass
                mmc_object.resonance.set_from(mmc_resonance)
                if mmc_mass > 0:
                    mmc_object.resonance_pt = mmc_resonance.Pt()
                mmc_object.MET = mmc_met.Mod()
                mmc_object.MET_x = mmc_met.X()
                mmc_object.MET_y = mmc_met.Y()
                mmc_object.MET_phi = math.pi - mmc_met.Phi()
                mmc_object.MET_vec.set_from(mmc_met)

            ############################
            # collinear and visible mass
            ############################
            vis_mass, collin_mass, tau1_x, tau2_x = mass.collinearmass(
                    tau1, tau2, METx, METy)

            tree.mass_vis_tau1_tau2 = vis_mass
            tree.mass_collinear_tau1_tau2 = collin_mass
            tau1.collinear_momentum_fraction = tau1_x
            tau2.collinear_momentum_fraction = tau2_x

            if datatype == datasets.MC and year == 2011:
                tree.ggf_weight = reweight_ggf(event, self.metadata.name)

            ###########################
            # Match jets to VBF partons
            ###########################
            if datatype == datasets.MC and 'VBF' in self.metadata.name and year == 2011:
                # get partons (already sorted by eta in hepmc) FIXME!!!
                parton1, parton2 = hepmc.get_VBF_partons(event)
                tree.mass_true_quark1_quark2 = (parton1.fourvect + parton2.fourvect).M()

                # order here needs to be revised since jets are no longer
                # sorted by eta but instead by pT
                PartonBlock.set(tree, parton1, parton2)
                if len(jets) >= 2:
                    jet1, jet2 = jets[:2]
                    for i, jet in zip((1, 2), (jet1, jet2)):
                        for parton in (parton1, parton2):
                            if utils.dR(jet.eta, jet.phi, parton.eta, parton.phi) < .8:
                                setattr(tree, 'jet%i_matched' % i, True)

            ###########################
            # truth matching
            ###########################
            if datatype == datasets.MC:
                # match only with visible true taus
                event.truetaus.select(
                        lambda tau: tau.vis_Et > 10 * GeV and abs(tau.vis_eta) < 2.5)

                if len(event.truetaus) > 2:
                    log.warning("too many true taus: %i" % len(event.truetaus))
                    for truetau in event.truetaus:
                        print "truth (pT: %.4f, eta: %.4f, phi: %.4f)" % (
                                truetau.pt, truetau.eta, truetau.phi),
                        if truetau.tauAssoc_index >= 0:
                            matched_tau = event.taus.getitem(truetau.tauAssoc_index)
                            print " ==> reco (pT: %.4f, eta: %.4f, phi: %.4f)" % (
                                    matched_tau.pt, matched_tau.eta, matched_tau.phi),
                            print "dR = %.4f" % truetau.tauAssoc_dr
                        else:
                            print ""
                    tree.error = True

                unmatched_reco = range(2)
                unmatched_truth = range(event.truetaus.len())
                matched_truth = []
                for i, tau in enumerate((tau1, tau2)):
                    matching_truth_index = tau.trueTauAssoc_index
                    if matching_truth_index >= 0:
                        unmatched_reco.remove(i)
                        # check that this tau / true tau was not previously matched
                        if (matching_truth_index not in unmatched_truth or
                            matching_truth_index in matched_truth):
                            log.warning("match collision!")
                            tau1.matched_collision = True
                            tau2.matched_collision = True
                            tree.truetau1_matched_collision = True
                            tree.truetau2_matched_collision = True
                            tree.error = True
                        else:
                            unmatched_truth.remove(matching_truth_index)
                            matched_truth.append(matching_truth_index)
                            tau.matched = True
                            tau.matched_dR = tau.trueTauAssoc_dr
                            setattr(tree, "truetau%i_matched" % (i+1), 1)
                            setattr(tree, "truetau%i_matched_dR" % (i+1),
                                    event.truetaus.getitem(
                                        matching_truth_index).tauAssoc_dr)
                            TrueTauBlock.set(tree, i+1,
                                    event.truetaus.getitem(matching_truth_index))

                for i, j in zip(unmatched_reco, unmatched_truth):
                    TrueTauBlock.set(tree, i+1, event.truetaus.getitem(j))

                tree.mass_vis_true_tau1_tau2 = (
                        tree.truetau1_fourvect_vis +
                        tree.truetau2_fourvect_vis).M()

            ###########################
            # Fill tau block
            ###########################

            # This must come after the RecoJetBlock is filled since
            # that sets the jet_beta for boosting the taus
            RecoTauBlock.set(event, tree, tau1, tau2)

            # TODO UPDATE:
            if validate:
                if datatype == datasets.MC:
                    print >> validate_log, event.mc_channel_number,
                print >> validate_log, event.RunNumber, event.EventNumber,
                print >> validate_log, "%.4f" % tree.pileup_weight,
                for idx in selected_idx:
                    print >> validate_log, idx, tree.tau_trigger_match_thresh[idx],
                print >> validate_log

                print "/" * 60
                print "/" * 60

                print "entry:", event._entry.value
                print "EventNumber:", event.EventNumber
                print
                print "Pileup weight:", tree.pileup_weight
                print
                print "trigger scale factors (taus ordered by decreasing pT):"
                for i, tau in enumerate(event.taus):
                    print
                    print "tau %d:" % (i + 1)
                    print "BDT ID Loose %d, Medium %d, Tight %d" % (
                        tau.JetBDTSigLoose,
                        tau.JetBDTSigMedium,
                        tau.JetBDTSigTight)
                    print "pT: %f" % tau.pt
                    print "matched trigger threshold: %d" % tau.trigger_match_thresh
                    print "matched trigger index: %d" % tau.trigger_match_index

                    loose = tau.trigger_eff_sf_loose
                    loose_high = tau.trigger_eff_sf_loose_high
                    loose_low = tau.trigger_eff_sf_loose_low

                    medium = tau.trigger_eff_sf_medium
                    medium_high = tau.trigger_eff_sf_medium_high
                    medium_low = tau.trigger_eff_sf_medium_low

                    tight = tau.trigger_eff_sf_tight
                    tight_high = tau.trigger_eff_sf_tight_high
                    tight_low = tau.trigger_eff_sf_tight_low

                    fmt = "%s: %f (high: %f low: %f)"
                    print fmt % ('loose', loose, loose_high, loose_low)
                    print fmt % ('medium', medium, medium_high, medium_low)
                    print fmt % ('tight', tight, tight_high, tight_low)
                print
                print "mass:"
                # print out values for comparing with Soshi's code
                vis_mass_alt, collin_mass_alt, tau1_x_alt, tau2_x_alt = \
                        mass.collinearmass_alt(tau1, tau2, METx, METy)
                print "vis_mass", "coll_mass", "x1", "x2"
                print vis_mass, collin_mass, tau1_x, tau2_x, "(me)"
                print vis_mass_alt, collin_mass_alt, tau1_x_alt, tau2_x_alt, "(Soshi)"
                print

            ######################
            # fill the output tree
            ######################
            tree.Fill(reset=True)

        self.output.cd()

        if validate:
            validate_log.close()
            # sort the output by event number for MC
            # sort +2 -3 -n skim2_validate_mc_125205.txt -o skim2_validate_mc_125205.txt

        ###############################################
        # flush any baskets remaining in memory to disk
        ###############################################
        tree.FlushBaskets()
        tree.Write()

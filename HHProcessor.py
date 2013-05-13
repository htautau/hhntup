import warnings
import numpy as np
warnings.filterwarnings('error', category=np.ComplexWarning)

# jet calibration sometimes gives jets with pT=0
#warnings.filterwarnings('error', 'transvers momentum = 0!', RuntimeWarning)

import ROOT
import math

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.math.physics.vector import Vector2, LorentzVector
from rootpy.plotting import Hist
from rootpy.extern.argparse import ArgumentParser

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.filtering import GRLFilter
from atlastools.batch import ATLASStudent

from higgstautau.mixins import *
from higgstautau import hepmc
from higgstautau import tautools
from higgstautau.models import *
from higgstautau.hadhad.models import *
from higgstautau import eventshapes
from higgstautau import eventview
from higgstautau.filters import *
from higgstautau.hadhad.filters import *
from higgstautau import mass
from higgstautau.mass import is_MET_bisecting
from higgstautau.embedding import EmbeddingPileupPatch, EmbeddingIsolation
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger.emulation import TauTriggerEmulation, update_trigger_trees
from higgstautau.trigger.matching import TauTriggerMatchIndex, TauTriggerMatchThreshold
from higgstautau.trigger.efficiency import TauTriggerEfficiency
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.overlap import TauJetOverlapRemoval
from higgstautau.patches import ElectronIDpatch, TauIDpatch
from higgstautau.pileup import (PileupReweight, PileupDataScale,
                                averageIntPerXingPatch,
                                get_pileup_reweighting_tool)
from higgstautau.hadhad.objects import define_objects
from higgstautau import log; log = log[__name__]

from goodruns import GRL
import subprocess

#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class HHProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def __init__(self, options, **kwargs):

        super(HHProcessor, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-terms', default=None)
        parser.add_argument('--redo-mmc', default=False, action='store_true')
        parser.add_argument('--student-verbose', default=False, action='store_true')
        self.args = parser.parse_args(options)
        if self.args.syst_terms is not None:
            self.args.syst_terms = set([
                eval('Systematics.%s' % term) for term in
                self.args.syst_terms.split(',')])

    @staticmethod
    def merge(inputs, output, metadata):

        # merge output trees
        root_output = output + '.root'
        subprocess.call(['hadd', root_output] + inputs)

        if metadata.datatype == datasets.DATA:
            # merge GRLs
            grl = GRL()
            for input in inputs:
                grl |= GRL('%s:/lumi' % input)
            grl.save('%s:/lumi' % root_output)

    def work(self):
        """
        This is the one function that all "ATLASStudent"s must implement.
        """
        datatype = self.metadata.datatype
        year = self.metadata.year
        redo_mmc = self.args.redo_mmc
        verbose = self.args.student_verbose

        # get pileup reweighting tool
        pileup_tool = get_pileup_reweighting_tool(
            year=year,
            use_defaults=True)

        OutputModel = RecoTauBlock + RecoJetBlock + EventVariables

        if datatype == datasets.MC:
            # only create truth branches for MC
            OutputModel += TrueTauBlock

            # add branches for VBF Higgs associated partons
            if 'VBF' in self.metadata.name:
                OutputModel += PartonBlock

        if datatype == datasets.EMBED:
            # add embedding systematics branches
            OutputModel += EmbeddingBlock

        onfilechange = []
        count_funcs = {}

        if datatype in (datasets.MC, datasets.EMBED):

            def mc_weight_count(event):
                return event.mc_event_weight

            count_funcs = {
                'mc_weight': mc_weight_count,
            }

        trigger_config = None

        if datatype != datasets.EMBED:
            # trigger config tool to read trigger info in the ntuples
            trigger_config = get_trigger_config()

            # update the trigger config maps on every file change
            onfilechange.append((update_trigger_config, (trigger_config,)))

        if datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file, tree):

                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))

        if datatype == datasets.DATA:
            merged_cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        else:
            merged_cutflow = Hist(2, 0, 2, name='cutflow', type='D')

        def update_cutflow(student, cutflow, name, file, tree):

            year = student.metadata.year
            datatype = student.metadata.datatype
            if datatype == datasets.MC:
                cutflow[0] += file.cutflow_event[0]
                cutflow[1] += file.cutflow_event_mc_weight[0]
            else:
                cutflow[0] += file.cutflow_event[0]

        onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                events=self.events,
                read_branches_on_demand=True,
                cache=True,
                onfilechange=onfilechange)

        # create output tree
        self.output.cd()
        tree = Tree(name='higgstautauhh', model=OutputModel)

        copied_variables = [
                'actualIntPerXing',
                'averageIntPerXing',
                'number_of_good_vertices',
                'RunNumber',
                'EventNumber',
                'lbn']

        tree.set_buffer(
                chain._buffer,
                branches=copied_variables,
                create_branches=True,
                visible=False)
        chain.always_read(copied_variables)

        # set the event filters
        event_filters = EventFilterList([
            GRLFilter(
                self.grl,
                passthrough=datatype not in (datasets.DATA, datasets.EMBED),
                count_funcs=count_funcs),
            # apply the PileupReweight before MCRunNumber
            PileupReweight(
                tool=pileup_tool,
                tree=tree,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            # set the RunNumber in MC to a random run number according to the
            # pileup tool
            MCRunNumber(
                pileup_tool=pileup_tool,
                passthrough=datatype != datasets.MC,
                count_funcs=count_funcs),
            #EmbeddingPileupPatch(
            #    passthrough=year > 2011 or datatype != datasets.EMBED,
            #    count_funcs=count_funcs),
            averageIntPerXingPatch(
                passthrough=year < 2012 or datatype != datasets.MC,
                count_funcs=count_funcs),
            #Triggers(
            #    year=year,
            #    old_skim=datatype == datasets.MC,
            #    passthrough=datatype == datasets.EMBED,
            #    count_funcs=count_funcs),
            #PriVertex(
            #    count_funcs=count_funcs),
            #LArError(
            #    count_funcs=count_funcs),
            # no need to recalibrate jets in 2012 (yet...)
            #JetCalibration(
            #    datatype=datatype,
            #    year=year,
            #    verbose=verbose,
            #    count_funcs=count_funcs),
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
            #LArHole(
            #    datatype=datatype,
            #    count_funcs=count_funcs),
            #JetCleaning(
            #    datatype=datatype,
            #    year=year,
            #    count_funcs=count_funcs),
            #ElectronVeto(
            #    count_funcs=count_funcs),
            #MuonVeto(
            #    year=year,
            #    count_funcs=count_funcs),
            #TauElectronVeto(2,
            #    count_funcs=count_funcs),
            #TauMuonVeto(2,
            #    count_funcs=count_funcs),
            #TauAuthor(2,
            #    count_funcs=count_funcs),
            #TauHasTrack(2,
            #    count_funcs=count_funcs),
            #TauPT(2,
            #    thresh=20 * GeV,
            #    count_funcs=count_funcs),
            #TauEta(2,
            #    count_funcs=count_funcs),
            #TauCrack(2,
            #    count_funcs=count_funcs),
            #TauLArHole(2,
            #    count_funcs=count_funcs),
            #TauTriggerMatchIndex(
            #    config=trigger_config,
            #    year=year,
            #    datatype=datatype,
            #    passthrough=datatype == datasets.EMBED,
            #    count_funcs=count_funcs),
            TauSelected(2,
                count_funcs=count_funcs),
            TauIDSelection(
                year=year,
                count_funcs=count_funcs),
            # apply this selection here since skim has lower threshold for data
            TauLeadSublead(
                lead=35 * GeV,
                sublead=25 * GeV,
                count_funcs=count_funcs),
            TruthMatching(
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            TauTriggerMatchThreshold(
                datatype=datatype,
                tree=tree,
                count_funcs=count_funcs),
            TauTriggerEfficiency(
                year=year,
                datatype=datatype,
                tree=tree,
                tes_systematic=self.args.syst_terms and (
                    Systematics.TES_TERMS & self.args.syst_terms),
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            PileupDataScale(
                year=year,
                passthrough=datatype not in (datasets.DATA, datasets.EMBED),
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
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            ggFReweighting(
                dsname=self.metadata.name,
                tree=tree,
                # no ggf reweighting for 2012 MC
                passthrough=datatype != datasets.MC or year != 2011,
                count_funcs=count_funcs),
            #TauTrackRecounting(
            #    year=year,
            #    count_funcs=count_funcs),
            MCWeight(
                datatype=datatype,
                tree=tree,
                passthrough=datatype == datasets.DATA,
                count_funcs=count_funcs),
            EmbeddingIsolation(
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

        self.filters['event'] = event_filters

        chain._filters += event_filters

        define_objects(chain, year, skim=False)

        # define tree objects
        tree.define_object(name='tau1', prefix='tau1_')
        tree.define_object(name='tau2', prefix='tau2_')
        tree.define_object(name='jet1', prefix='jet1_')
        tree.define_object(name='jet2', prefix='jet2_')

        """ Associations not currently implemented in rootpy
        chain.define_association(origin='taus', target='truetaus', prefix='trueTauAssoc_', link='index')
        chain.define_association(origin='truetaus', target='taus', prefix='tauAssoc_', link='index')
        """

        if redo_mmc:
            log.info("will recalculate MMC output")
            # create MMC object
            mmc = mass.MMC(year=year, channel='hh')

        # entering the main event loop...
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

            ##########################
            # Jet and sum pt variables
            ##########################
            tree.numJets = len(event.jets)
            # sum pt with only the two leading jets
            tree.sum_pt = sum(
                    [tau1.pt, tau2.pt] +
                    [jet.pt for jet in jets[:2]])
            # sum pt with all selected jets
            tree.sum_pt_full = sum(
                    [tau1.pt, tau2.pt] +
                    [jet.pt for jet in jets])
            # vector sum pt with two leading jets and MET
            tree.vector_sum_pt = sum(
                    [tau1.fourvect, tau2.fourvect] +
                    [jet.fourvect for jet in jets[:2]] +
                    [MET_4vect]).Pt()

            ##########################
            # MMC Mass
            ##########################
            if redo_mmc:
                mmc_mass, mmc_resonance, mmc_met = mmc.mass(
                        tau1, tau2,
                        METx, METy, sumET,
                        len(event.jets),
                        method=0)
                if verbose:
                    log.info("MMC: %f" % mmc_mass)
            else:
                # use MMC values from skim
                mmc_mass = event.tau_MMC_mass
                mmc_resonance = event.tau_MMC_resonance
                mmc_met = Vector2(event.tau_MMC_MET_x, event.tau_MMC_MET_y)

            tree.mass_mmc_tau1_tau2 = mmc_mass
            tree.mmc_resonance.set_from(mmc_resonance)
            if mmc_mass > 0:
                tree.mmc_resonance_pt = mmc_resonance.Pt()
            tree.MET_mmc = mmc_met.Mod()
            tree.MET_mmc_x = mmc_met.X()
            tree.MET_mmc_y = mmc_met.Y()
            tree.MET_mmc_phi = math.pi - mmc_met.Phi()
            tree.MET_mmc_vec.set_from(mmc_met)

            mass_vis, mass_col, tau1_x, tau2_x = mass.collinearmass(
                    tau1, tau2, METx, METy)
            tree.mass_collinear_tau1_tau2 = mass_col
            tree.mass_vis_tau1_tau2 = mass_vis
            tree.tau1_x = tau1_x
            tree.tau2_x = tau2_x
            tree.tau_x_product = tau1_x * tau2_x
            tree.tau_x_sum = tau1_x + tau2_x

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
                            tree.trueTau1_matched_collision = True
                            tree.trueTau2_matched_collision = True
                            tree.error = True
                        else:
                            unmatched_truth.remove(matching_truth_index)
                            matched_truth.append(matching_truth_index)
                            tau.matched = True
                            tau.matched_dR = tau.trueTauAssoc_dr
                            setattr(tree, "trueTau%i_matched" % (i+1), 1)
                            setattr(tree, "trueTau%i_matched_dR" % (i+1),
                                    event.truetaus.getitem(
                                        matching_truth_index).tauAssoc_dr)
                            TrueTauBlock.set(tree, i+1,
                                    event.truetaus.getitem(matching_truth_index))

                for i, j in zip(unmatched_reco, unmatched_truth):
                    TrueTauBlock.set(tree, i+1, event.truetaus.getitem(j))

                tree.mass_vis_true_tau1_tau2 = (
                        tree.trueTau1_fourvect_vis +
                        tree.trueTau2_fourvect_vis).M()

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

            # Fill tau block
            # This must come after the RecoJetBlock is filled since
            # that sets the jet_beta for boosting the taus
            RecoTauBlock.set(event, tree, tau1, tau2)

            # Fill output ntuple
            tree.Fill(reset=True)

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

        if datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        merged_cutflow.Write()

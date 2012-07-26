import ROOT
import math

from argparse import ArgumentParser

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.tree.cutflow import Cutflow
from rootpy.math.physics.vector import Vector2
from rootpy.plotting import Hist
from rootpy.io import open as ropen

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
from higgstautau.hadhad.categories import *
from higgstautau import mass
#from higgstautau.mass.ditaumass import HAD1P, HAD3P
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.pileup import PileupReweighting, TPileupReweighting
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration
from higgstautau.overlap import TauJetOverlapRemoval

from goodruns import GRL
import subprocess

from externaltools import TauFakeRates
from ROOT import TauFakeRates as TFR

#ROOT.gErrorIgnoreLevel = ROOT.kFatal
YEAR = 2011
VERBOSE = False

class HHProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def __init__(self, options, **kwargs):

        super(HHProcessor, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-terms', default=None)
        self.args = parser.parse_args(options)
        if self.args.syst_terms is not None:
            self.args.syst_terms = [
                eval(term) for term in
                self.args.syst_terms.split(',')]

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

        # fake rate scale factor tool
        fakerate_table = TauFakeRates.get_resource('FakeRateScaleFactor.txt')
        fakerate_tool = TFR.FakeRateScaler(fakerate_table)

        # trigger config tool to read trigger info in the ntuples
        trigger_config = get_trigger_config()

        OutputModel = RecoTauBlock + RecoJetBlock + EventVariables

        if self.metadata.datatype == datasets.MC:
            # only create truth branches for MC
            OutputModel += TrueTauBlock

            # add branches for VBF Higgs associated partons
            if 'VBF' in self.metadata.name:
                OutputModel += PartonBlock

        onfilechange = []
        if self.metadata.datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file, tree):

                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))

        # update the trigger config maps on every file change
        onfilechange.append((update_trigger_config, (trigger_config,)))

        if self.metadata.datatype == datasets.DATA:
            merged_cutflow = Hist(1, 0, 1, name='cutflow', type='D')
        else:
            merged_cutflow = Hist(2, 0, 2, name='cutflow', type='D')

        def update_cutflow(student, cutflow, name, file, tree):

            cutflow += file.cutflow

        onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        chain = TreeChain(self.metadata.treename,
                         files=self.files,
                         events=self.events,
                         cache=True,
                         cache_size=10000000,
                         learn_entries=30,
                         onfilechange=onfilechange)

        # create output tree
        self.output.cd()
        tree = Tree(name='higgstautauhh', model=OutputModel)

        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']

        tree.set_buffer(
                chain.buffer,
                branches=copied_variables,
                create_branches=True,
                visible=False)
        chain.always_read(copied_variables)

        # set the event filters
        event_filters = EventFilterList([
            GRLFilter(
                self.grl,
                passthrough=self.metadata.datatype != datasets.DATA),
            Triggers(
                datatype=self.metadata.datatype,
                year=YEAR,
                skim=False),
            JetCalibration(
                year=YEAR,
                datatype=self.metadata.datatype,
                verbose=VERBOSE),
            # PUT THE SYSTEMATICS "FILTER" BEFORE
            # ANY FILTERS THAT REFER TO OBJECTS
            # BUT AFTER CALIBRATIONS
            Systematics(
                terms=self.args.syst_terms,
                year=YEAR,
                datatype=self.metadata.datatype,
                verbose=VERBOSE),
            # since the jet recalibration is applied the MET must be
            # recalculated even if no other systematics are applied.
            PriVertex(),
            LArError(),
            LArHole(datatype=self.metadata.datatype),
            JetCleaning(),
            #JetCrackVeto(),
            ElectronVeto(),
            MuonVeto(),
            TauAuthor(2),
            TauHasTrack(2),
            TauMuonVeto(2),
            TauElectronVeto(2),
            TauPT(2),
            TauEta(2),
            TauCrack(2),
            TauLArHole(2), # only veto taus, not entire event
            TauIDMedium(2),
            TauTriggerMatch(
                config=trigger_config,
                year=YEAR,
                datatype=self.metadata.datatype,
                skim=False,
                tree=tree),
            TauLeadSublead(
                lead=35*GeV,
                sublead=25*GeV),
            JetSelection(),
            TauJetOverlapRemoval(),
        ])

        self.filters['event'] = event_filters

        chain.filters += event_filters

        cutflow = Cutflow()

        # define tree collections
        chain.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        chain.define_collection(name="taus_EF", prefix="trig_EF_tau_",
                                size="trig_EF_tau_n", mix=TauFourMomentum)

        # jet_* etc. is AntiKt4LCTopo_* in tau-perf D3PDs
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        chain.define_collection(name="truetaus", prefix="trueTau_", size="trueTau_n", mix=MCTauFourMomentum)
        chain.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        chain.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n")
        chain.define_collection(name="electrons", prefix="el_", size="el_n")
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        # define tree objects
        tree.define_object(name='tau1', prefix='tau1_')
        tree.define_object(name='tau2', prefix='tau2_')
        tree.define_object(name='jet1', prefix='jet1_')
        tree.define_object(name='jet2', prefix='jet2_')

        """ Associations not currently implemented in rootpy
        chain.define_association(origin='taus', target='truetaus', prefix='trueTauAssoc_', link='index')
        chain.define_association(origin='truetaus', target='taus', prefix='tauAssoc_', link='index')
        """

        if self.metadata.datatype == datasets.MC:
            # Initialize the pileup reweighting tool
            pileup_tool = TPileupReweighting()
            #pileup_tool.AddConfigFile('/global/endw/mc11_7TeV/higgs_tautau_hh_reskim_p851/TPileupReweighting.prw.root')
            pileup_tool.AddConfigFile('higgstautau/pileup/mc11c_defaults.prw.root')
            pileup_tool.AddLumiCalcFile('grl/2011/lumicalc/hadhad/ilumicalc_histograms_None_178044-191933.root')
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetUnrepresentedDataAction(1)
            pileup_tool.Initialize()

        # entering the main event loop...
        for event in chain:
            tree.reset()
            cutflow.reset()

            # taus are already sorted by pT in TauLeadSublead filter
            tau1, tau2 = event.taus

            jets = list(event.jets)
            # sort by decreasing pT
            jets.sort(key=lambda jet: jet.pt, reverse=True)
            leading_jets = []

            current_channel = CATEGORY_GGF
            # leading jet above 50 GeV
            if jets and jets[0].pt > 50 * GeV:
                leading_jets.append(jets[0])
                current_channel = CATEGORY_BOOSTED
                # subleading jet above 30
                if len(jets) >= 2 and jets[1].pt > 30 * GeV:
                    leading_jets.append(jets[1])
                    current_channel = CATEGORY_VBF
            tree.category = current_channel

            if current_channel == CATEGORY_VBF: # VBF optimized
                jet1, jet2 = leading_jets
                RecoJetBlock.set(tree, jet1, jet2)
                """
                Reco tau variables
                This must come after the RecoJetBlock is filled since
                that sets the jet_beta for boosting the taus
                """
                sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                        [tau1.fourvect,
                         tau2.fourvect,
                         jet1.fourvect,
                         jet2.fourvect])

                # sphericity
                tree.sphericity = sphericity
                # aplanarity
                tree.aplanarity = aplanarity

                sphericity_full, aplanarity_full = eventshapes.sphericity_aplanarity(
                        [tau1.fourvect,
                         tau2.fourvect] + [jet.fourvect for jet in jets])

                # boosted sphericity
                tree.sphericity_full = sphericity_full
                # boosted aplanarity
                tree.aplanarity_full = aplanarity_full

                sphericity_b, aplanarity_b = eventshapes.sphericity_aplanarity(
                        [tau1.fourvect_boosted,
                         tau2.fourvect_boosted,
                         jet1.fourvect_boosted,
                         jet2.fourvect_boosted])

                # boosted sphericity
                tree.sphericity_boosted = sphericity_b
                # boosted aplanarity
                tree.aplanarity_boosted = aplanarity_b

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

            elif current_channel == CATEGORY_BOOSTED:
                jet1 = leading_jets[0]
                RecoJetBlock.set(tree, jet1)
                """
                Reco tau variables
                This must come after the RecoJetBlock is filled since
                that sets the jet_beta for boosting the taus
                """
                sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                        [tau1.fourvect,
                         tau2.fourvect,
                         jet1.fourvect])

                # sphericity
                tree.sphericity = sphericity
                # aplanarity
                tree.aplanarity = aplanarity

                sphericity_full, aplanarity_full = eventshapes.sphericity_aplanarity(
                        [tau1.fourvect,
                         tau2.fourvect] + [jet.fourvect for jet in jets])

                # boosted sphericity
                tree.sphericity_full = sphericity_full
                # boosted aplanarity
                tree.aplanarity_full = aplanarity_full

            # Jet variables
            tree.numJets = len(event.jets)
            tree.sum_pt = sum([tau1.pt, tau2.pt] +
                              [jet.pt for jet in leading_jets])
            tree.sum_pt_full = sum([tau1.pt, tau2.pt] +
                                   [jet.pt for jet in jets])

            # MET
            METx = event.MET_RefFinal_BDTMedium_etx
            METy = event.MET_RefFinal_BDTMedium_ety
            MET_vect = Vector2(METx, METy)
            MET_3vect = Vector3(METx, METy, 0.)
            MET = event.MET_RefFinal_BDTMedium_et
            tree.MET = MET
            tree.MET_phi = event.MET_RefFinal_BDTMedium_phi
            sumET = event.MET_RefFinal_BDTMedium_sumet
            tree.HT = sumET
            tree.MET_sig = (2. * MET / GeV) / (utils.sign(sumET) * sqrt(abs(sumET / GeV)))
            MET_res = 6.14 * math.sqrt(GeV) + 0.5 * math.sqrt(abs(sumET))

            tau1_2vector = Vector2(tau1.fourvect.Px(), tau1.fourvect.Py())
            tau2_2vector = Vector2(tau2.fourvect.Px(), tau2.fourvect.Py())
            tree.MET_centrality = eventshapes.phi_centrality(tau1_2vector,
                                                             tau2_2vector,
                                                             MET_vect)

            # Mass
            mmc_mass, mmc_pt, mmc_met = mass.missingmass(tau1, tau2, METx, METy, sumET)
            tree.mass_mmc_tau1_tau2 = mmc_mass
            tree.higgs_pt = mmc_pt
            tree.MET_mmc = mmc_met

            """
            if tau1.numTrack <= 1:
                taumode1 = HAD1P
            else:
                taumode1 = HAD3P

            if tau2.numTrack <= 1:
                taumode2 = HAD1P
            else:
                taumode2 = HAD3P

            tree.mass_dtm_tau1_tau2 = mass.ditaumass(tau1.fourvect, taumode1,
                                                     tau2.fourvect, taumode2,
                                                     METx, METy, MET_res) / GeV
            tree.mass_dtm_tau1_tau2_scan = mass.ditaumass_scan(tau1.fourvect, taumode1,
                                                     tau2.fourvect, taumode2,
                                                     METx, METy, MET_res, 5) / GeV
            """
            collin_mass, tau1_x, tau2_x = mass.collinearmass(tau1, tau2, METx, METy)
            tree.mass_collinear_tau1_tau2 = collin_mass
            tree.tau1_x = tau1_x
            tree.tau2_x = tau2_x

            tree.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or
                                    (vtx.type == 3 and vtx.nTracks >= 2)])

            # Match jets to VBF partons
            if self.metadata.datatype == datasets.MC:
                if 'VBF' in self.metadata.name:
                    # get partons (already sorted by eta in hepmc)
                    parton1, parton2 = hepmc.get_VBF_partons(event)
                    tree.mass_true_quark1_quark2 = (parton1.fourvect + parton2.fourvect).M()

                    # order here needs to be revised since jets are no longer
                    # sorted by eta but instead by pT
                    PartonBlock.set(tree, parton1, parton2)
                    if current_channel == CATEGORY_VBF:
                        for i, jet in zip((1, 2), (jet1, jet2)):
                            for parton in (parton1, parton2):
                                if utils.dR(jet.eta, jet.phi, parton.eta, parton.phi) < .8:
                                    setattr(tree, 'jet%i_matched' % i, True)

            # Truth-matching
            if self.metadata.datatype == datasets.MC:
                # match only with visible true taus
                event.truetaus.select(lambda tau: tau.vis_Et > 10 * GeV and abs(tau.vis_eta) < 2.5)

                if len(event.truetaus) > 2:
                    print "ERROR: too many true taus: %i" % len(event.truetaus)
                    for truetau in event.truetaus:
                        print "truth (pT: %.4f, eta: %.4f, phi: %.4f)" % (truetau.pt, truetau.eta, truetau.phi),
                        if truetau.tauAssoc_index >= 0:
                            matched_tau = event.taus.getitem(truetau.tauAssoc_index)
                            print " ==> reco (pT: %.4f, eta: %.4f, phi: %.4f)" % (matched_tau.pt, matched_tau.eta, matched_tau.phi),
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
                        if matching_truth_index not in unmatched_truth or \
                           matching_truth_index in matched_truth:
                            print "ERROR: match collision!"
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
                            setattr(tree, "trueTau%i_matched_dR" % (i+1), event.truetaus.getitem(matching_truth_index).tauAssoc_dr)
                            TrueTauBlock.set(tree, i+1, event.truetaus.getitem(matching_truth_index))

                for i, j in zip(unmatched_reco, unmatched_truth):
                    TrueTauBlock.set(tree, i+1, event.truetaus.getitem(j))

                tree.mass_vis_true_tau1_tau2 = (tree.trueTau1_fourvect_vis + tree.trueTau2_fourvect_vis).M()

                for tau in (tau1, tau2):
                    if tau.matched:
                        # efficiency scale factor
                        tau.weight = 1.
                    else:
                        # fake rate scale factor
                        if event.RunNumber >= 188902:
                            trig = "EF_tau%dT_medium1"
                        else:
                            trig = "EF_tau%d_medium1"
                        tau.weight = fakerate_tool.getScaleFactor(
                                tau.pt, "Medium",
                                trig % tau.trigger_match_thresh)

            # fill tau block
            RecoTauBlock.set(event, tree, tau1, tau2)

            # fill output ntuple
            tree.cutflow = cutflow.int()
            if self.metadata.datatype == datasets.MC:
                # set the event weight
                tree.pileup_weight = pileup_tool.GetCombinedWeight(event.RunNumber,
                                                                   event.mc_channel_number,
                                                                   event.averageIntPerXing)
                tree.mc_weight = event.mc_event_weight
            tree.Fill()

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

        if self.metadata.datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        merged_cutflow.Write()

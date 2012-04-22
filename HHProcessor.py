import ROOT
import math

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
from higgstautau.hadhad.extrafilters import *
from higgstautau import mass
#from higgstautau.mass.ditaumass import HAD1P, HAD3P
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.pileup import PileupReweighting, TPileupReweighting

from goodruns import GRL
import subprocess


#ROOT.gErrorIgnoreLevel = ROOT.kFatal


class HHProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

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

        merged_cutflow = Hist(1, 0, 1, name='cutflow', type='D')
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
        tree_hh_2jet = Tree(name=self.metadata.name + '_2jet', model=OutputModel)
        #tree_hh_1jet = Tree(name=self.metadata.name + '_1jet', model=OutputModel)
        tree_hh_01jet = Tree(name=self.metadata.name + '_01jet', model=OutputModel)

        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']

        if self.metadata.datatype == datasets.MC:
            copied_variables += mc_triggers

        tree_hh_2jet.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)
        tree_hh_01jet.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)

        #tree_hh_1jet.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)
        #tree_hh_0jet.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)

        chain.always_read(copied_variables)

        # set the event filters
        # passthrough for MC for trigger acceptance studies
        event_filters = EventFilterList([
            GRLFilter(self.grl, passthrough=self.metadata.datatype != datasets.DATA),
            Triggers(),
            PriVertex(),
            LArError(),
            #SomeJets(),
            #SomeTaus(),
            LArHole(datatype=self.metadata.datatype),
            JetCleaning(),
            #JetCrackVeto(),
            ElectronVeto(),
            MuonVeto(),
            TauAuthor(),
            TauHasTrack(),
            TauCharge(),
            TauMuonVeto(),
            TauElectronVeto(),
            TauPT(),
            TauEta(),
            TauCrack(),
            #TauLArHole(), #only veto taus, not entire event
            TauLoose(),
            TauTriggerMatch(config=trigger_config),
            TauLeadSublead(lead=35*GeV,
                           sublead=25*GeV),
            #TauJVF(),
            #Tau1Track3Track()
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
        tree_hh_2jet.define_object(name='tau1', prefix='tau1_')
        tree_hh_2jet.define_object(name='tau2', prefix='tau2_')
        tree_hh_2jet.define_object(name='jet1', prefix='jet1_')
        tree_hh_2jet.define_object(name='jet2', prefix='jet2_')

        tree_hh_01jet.define_object(name='tau1', prefix='tau1_')
        tree_hh_01jet.define_object(name='tau2', prefix='tau2_')
        #tree_hh_1jet.define_object(name='jet1', prefix='jet1_')

        #tree_hh_0jet.define_object(name='tau1', prefix='tau1_')
        #tree_hh_0jet.define_object(name='tau2', prefix='tau2_')

        """ Associations not currently implemented in rootpy
        chain.define_association(origin='taus', target='truetaus', prefix='trueTauAssoc_', link='index')
        chain.define_association(origin='truetaus', target='taus', prefix='tauAssoc_', link='index')
        """

        #CHAN_0JET, CHAN_1JET, CHAN_2JET = range(3)
        CHAN_01JET, CHAN_2JET = range(2)

        if self.metadata.datatype == datasets.MC:
            # Initialize the pileup reweighting tool
            pileup_tool = TPileupReweighting()
            pileup_tool.AddConfigFile(PileupReweighting.get_resource('%s_defaults.prw.root' % self.metadata.category))
            pileup_tool.AddLumiCalcFile('grl/lumicalc/hadhad/ilumicalc_histograms_None_178044-191933.root')
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetUnrepresentedDataAction(1)
            pileup_tool.Initialize()

        # entering the main event loop...
        for event in chain:

            #tree_hh_0jet.reset()
            #tree_hh_1jet.reset()
            tree_hh_01jet.reset()
            tree_hh_2jet.reset()

            cutflow.reset()

            # remove pileup and UE
            # event.jets.select(lambda jet: True if abs(jet.eta) > 2.1 else jet.jvtxf > .5)

            # find pair of jets with highest dijet mass
            """
            highest_mass = 0
            best_jets = ()
            for i, jet1 in enumerate(event.jets[:-1]):
                for jet2 in event.jets[i+1:]:
                    dijet_mass = (jet1.fourvect + jet2.fourvect).M()
                    if dijet_mass > highest_mass:
                        highest_mass = dijet_mass
                        best_jets = (jet1, jet2)
            """
            # Sort the taus by BDT score
            # event.taus.sort(key=lambda tau: tau.BDTJetScore, reverse=True)
            # Take the two taus with the highest BDT score
            # taus = event.taus[:2]
            tau1, tau2 = event.taus
            # Jet selection
            event.jets.select(lambda jet: jet.pt > 25 * GeV and abs(jet.eta) < 4.5)

            # remove overlap with taus
            event.jets.select(lambda jet: not any([tau for tau in event.taus if \
                                                   (utils.dR(jet.eta, jet.phi, tau.eta, tau.phi) < .2)]))

            # select VBF jets
            jets = list(event.jets)
            # sort by decreasing pT
            jets.sort(key=lambda jet: jet.pt, reverse=True)
            # take the two jets with the highest pT
            leading_jets = jets[:2]
            # sort by increasing eta
            leading_jets.sort(key=lambda jet: jet.eta)

            """ VBF cuts proposed by Zinonas
            if len(leading_jets) >= 2:
                jet1, jet2 = leading_jets
                # pT of leading and subleading jets
                if jet1.pt > 40 * GeV and jet2.pt > 30 * GeV:
                    # opposite hemispheres
                    if jet1.eta * jet2.eta < 0:
                        # eta gap requirement
                        if abs(jet1.eta - jet2.eta) > 3:
                            # mass requirement
                            if (jet1.fourvect + jet2.fourvect).M() > 200 * GeV:
                                # require that taus are in the middle
                                if jet1.eta < tau1.eta < jet2.eta and \
                                   jet1.eta < tau2.eta < jet2.eta:
                                       VBF_jets = [jet1, jet2]
            """

            if len(leading_jets) >= 2: # VBF optimized
                current_channel = CHAN_2JET
                tree = tree_hh_2jet
                jet1, jet2 = leading_jets
                RecoJetBlock.set(tree, jet1, jet2)
                """
                Reco tau variables
                This must come after the RecoJetBlock is filled since
                that sets the jet_beta for boosting the taus
                """
                RecoTauBlock.set(event, tree, tau1, tau2)
                sphericity, aplanarity = eventshapes.sphericity_aplanarity(
                        [tau1.fourvect,
                         tau2.fourvect,
                         jet1.fourvect,
                         jet2.fourvect])

                # sphericity
                tree.sphericity = sphericity
                # aplanarity
                tree.aplanarity = aplanarity

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
                tree.tau1_centrality = eventshapes.eta_centrality(tau1.fourvect.Eta(),
                                                                  jet1.fourvect.Eta(),
                                                                  jet2.fourvect.Eta())

                tree.tau2_centrality = eventshapes.eta_centrality(tau2.fourvect.Eta(),
                                                                  jet1.fourvect.Eta(),
                                                                  jet2.fourvect.Eta())
                # boosted tau centrality
                tree.tau1_centrality_boosted = eventshapes.eta_centrality(
                        tau1.fourvect_boosted.Eta(),
                        jet1.fourvect_boosted.Eta(),
                        jet2.fourvect_boosted.Eta())

                tree.tau2_centrality_boosted = eventshapes.eta_centrality(
                        tau2.fourvect_boosted.Eta(),
                        jet1.fourvect_boosted.Eta(),
                        jet2.fourvect_boosted.Eta())

                """
                elif len(leading_jets) == 1: # one jet
                current_channel = CHAN_1JET
                tree = tree_hh_1jet
                RecoTauBlock.set(event, tree, tau1, tau2)
                """

            else: # 0 or 1 jet
                current_channel = CHAN_01JET
                tree = tree_hh_01jet
                RecoTauBlock.set(event, tree, tau1, tau2)


            """
            Jet variables
            """
            tree.numJets = len(event.jets)

            """
            MET
            """
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


            """
            Mass
            """
            tree.mass_mmc_tau1_tau2 = mass.missingmass(tau1, tau2, METx, METy, sumET)

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

            """
            Experimenting here....
            match jets to VBF jets
            """
            if self.metadata.datatype == datasets.MC:
                if 'VBF' in self.metadata.name:
                    # get partons (already sorted by eta in hepmc)
                    parton1, parton2 = hepmc.get_VBF_partons(event)
                    tree.mass_true_quark1_quark2 = (parton1.fourvect + parton2.fourvect).M()
                    PartonBlock.set(tree, parton1, parton2)
                    if current_channel == CHAN_2JET:
                        for i, jet in zip((1, 2), (jet1, jet2)):
                            for parton in (parton1, parton2):
                                if utils.dR(jet.eta, jet.phi, parton.eta, parton.phi) < .8:
                                    setattr(tree, 'jet%i_matched' % i, True)

            """
            Truth-matching
            presently not possible in SMWZ D3PDs
            """
            if self.metadata.datatype == datasets.MC and hasattr(event, "trueTau_n"):
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
                            tree.tau1_matched_collision = True
                            tree.tau2_matched_collision = True
                            tree.trueTau1_matched_collision = True
                            tree.trueTau2_matched_collision = True
                            tree.error = True
                        else:
                            unmatched_truth.remove(matching_truth_index)
                            matched_truth.append(matching_truth_index)
                            setattr(tree, "tau%i_matched" % (i+1), 1)
                            setattr(tree, "tau%i_matched_dR" % (i+1), tau.trueTauAssoc_dr)
                            setattr(tree, "trueTau%i_matched" % (i+1), 1)
                            setattr(tree, "trueTau%i_matched_dR" % (i+1), event.truetaus.getitem(matching_truth_index).tauAssoc_dr)
                            TrueTauBlock.set(tree, i+1, event.truetaus.getitem(matching_truth_index))

                for i, j in zip(unmatched_reco, unmatched_truth):
                    TrueTauBlock.set(tree, i+1, event.truetaus.getitem(j))

                tree.mass_vis_true_tau1_tau2 = (tree.trueTau1_fourvect_vis + tree.trueTau2_fourvect_vis).M()

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
        tree_hh_2jet.FlushBaskets()
        tree_hh_2jet.Write()
        #tree_hh_1jet.FlushBaskets()
        #tree_hh_1jet.Write()
        tree_hh_01jet.FlushBaskets()
        tree_hh_01jet.Write()

        if self.metadata.datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        merged_cutflow.Write()

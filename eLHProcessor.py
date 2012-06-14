import ROOT
import math

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.tree.cutflow import Cutflow
from rootpy.math.physics.vector import Vector2, LorentzVector
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
from higgstautau.lephad.models import *
from higgstautau import eventshapes
from higgstautau import eventview
from higgstautau.filters import *
from higgstautau.lephad.filters import *
from higgstautau.lephad.correctiontools import *
#from higgstautau.filters import PriVertex, JetCleaning, LArError, LArHole
from higgstautau import mass
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger import utils as triggerutils
from higgstautau.pileup import TPileupReweighting

from goodruns import GRL
import subprocess

import random

class eLHProcessor(ATLASStudent):
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

        OutputModel = RecoTauElectronBlock + TauElectronEventVariables + RecoMET

        onfilechange = []
        if self.metadata.datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file, tree):
                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))

        # update the trigger config maps on every file change
        onfilechange.append((update_trigger_config, (trigger_config,)))

        merged_cutflow = Hist(7, 0, 7, name='cutflow', type='D')

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
        tree_train = Tree(name=self.metadata.name + '_elh_train', model=OutputModel)
        tree_test = Tree(name=self.metadata.name + '_elh_test', model=OutputModel)

        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']

        tree_train.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)
        tree_test.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)

        chain.always_read(copied_variables)

        ###########################
        ## Set the event filters ##
        ###########################

        # Use single-lepton triggers
        Trigger = eMCSLTriggers
        if ( self.metadata.datatype == datasets.DATA ):
            Trigger = eSLTriggers

        # passthrough for MC for trigger acceptance studies
        event_filters = EventFilterList([
                SetElectronsFourVector(),
                # Require the event to fire a single-lepton trigger.
                Trigger(),
                GRLFilter(self.grl, passthrough=self.metadata.datatype != datasets.DATA),
                # At least 1 PV (not pile-up vertex) with Ntracks > 3.
                PriVertex(),
                MuonPtSmearing(datatype=self.metadata.datatype),
                EgammaERescaling(datatype=self.metadata.datatype),
                # Keep only jets with pt > 20 GeV. Keep event.
                JetPreSelection(),
                # Keep only muons with pt > 10 GeV, |eta| < 2.5, loose ID, has_good_track. Keep event.
                MuonPreSelection(),
                # Keep only electrons with cl_Et > 15 GeV, eta acceptance, medium ID, author. Keep event.
                ElectronPreSelection(),
                # Remove jets matching muons or electrons.
                JetOverlapRemoval(),
                # Remove bad jets.
                JetCleaning(eta_max = 9999.0),
                # Event should not pass if there is an electron in the LArHole region in the corresponding run #'s.
                ElectronLArHole(),
                # Event should not pass if there is a tau in the LArHole region in the corresponding run #'s.
                TauLArHole(),
                # Event should not pass if there is a jet (with enough high pt) in the LArHole region in the corresponding run #'s.
                LArHole(datatype=self.metadata.datatype),
                # Event should not pass if larError flag is > 1.
                LArError(),
                # Remove electrons matching muons.
                LeptonOverlapRemoval(),
                # There should be only 1 electron or muon.
                DileptonVeto(),
                # Remove event if there is not at least one tightPP electron with cl_Et > 25 GeV.
                ElectronSelection(),
                # There must be at least one tau with pt > 20 GeV, |eta| < 2.5, numTrack = 1 or 3, BDT medium ID, charge = +/-1, author != 2, EleBDTMedium == 0, muonVeto == 0.
                TauPreSelection(),
                # There must be only one tau.
                TauSelection(),
                # Keep only jets with pt > 25 GeV, |eta| < 4.5, JVF > 0.75 if |eta| < 2.4.
                JetSelection(),
                # Remove taus that overlap with selected muons or electrons, and remove jets that overlap with surviving selected taus.
                FinalOverlapRemoval()
                ])

        self.filters['event'] = event_filters

        chain.filters += event_filters

        cutflow = Cutflow()

        # define tree collections
        chain.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n", mix=FourMomentum)
        chain.define_collection(name="electrons", prefix="el_", size="el_n", mix=ElFourMomentum)
        chain.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        chain.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        # define tree objects
        tree_train.define_object(name='tau', prefix='tau_')
        tree_train.define_object(name='electron', prefix='el_')
        tree_test.define_object(name='tau', prefix='tau_')
        tree_test.define_object(name='electron', prefix='el_')

        if self.metadata.datatype == datasets.MC:
            # Initialize the pileup reweighting tool
            pileup_tool = TPileupReweighting()
            pileup_tool.AddConfigFile('higgstautau/pileup/%s_defaults.prw.root' % self.metadata.category)   # <--- Default file, NOT the one generated with the tool
            #pileup_tool.AddLumiCalcFile('higgstautau/lephad/external/Pileup/ilumicalc_histograms_None_178044-191933_slimmed.root')
            pileup_tool.AddLumiCalcFile('grl/2011/lumicalc/lephad/ilumicalc_histograms_None_178044-191933.root') # <--- "None" ? Shouldn't be with a trigger ?
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetUnrepresentedDataAction(2)
            pileup_tool.Initialize()

        # entering the main event loop...
        for event in chain: # <--- the filters are applied at this point; the loop is over events that pass the filters

            tree_train.reset()
            tree_test.reset()
            cutflow.reset()

            # Select if the event goes into the training or the testing tree
            tree = None
            #if event.EventNumber % 2 == 0:
            #    tree = tree_train
            #else:
            #    tree = tree_test
            if random.random() < 0.5:
                tree = tree_train
            else:
                tree = tree_test

            # Select tau with highest BDT score and surviving electron
            Tau = event.taus[0]
            Electron = event.electrons[0]

            """
            RecoTauElectronBlock filling
            """
            RecoTauElectronBlock.set(event, tree, Tau, Electron)

            """
            Jets
            """
            numJets = len(event.jets)
            tree.numJets = numJets

            numJets30 = 0
            numJets35 = 0

            for jet in event.jets:
                tree.jet_fourvect.push_back(jet.fourvect)
                tree.jet_jvtxf.push_back(jet.jvtxf)
                tree.jet_btag.push_back(jet.flavor_weight_JetFitterCOMBNN)
                if jet.fourvect.Pt() > 30*GeV:
                    numJets30 += 1
                    if jet.fourvect.Pt() > 30*GeV:
                        numJets35 += 1

            tree.numJets30 = numJets30
            tree.numJets35 = numJets35


            """
            Miscellaneous
            """
            tree.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or
                                         (vtx.type == 3 and vtx.nTracks >= 2)])


            # HT
            HT = 0
            HT += Tau.fourvect.Pt()
            HT += Electron.fourvect.Pt()
            for jet in event.jets:
                HT += jet.fourvect.Pt()
            tree.HT = HT


            # missing ET
            METx = event.MET_RefFinal_BDTMedium_etx
            METy = event.MET_RefFinal_BDTMedium_ety
            MET_vect = Vector2(METx, METy)
            MET = MET_vect.Mod()
            tree.MET = MET
            getattr(tree, 'MET_vect').set_from(MET_vect)
            sumET = event.MET_RefFinal_BDTMedium_sumet
            tree.MET_sig = (2. * MET_vect.Mod() / GeV) / (utils.sign(sumET) * sqrt(abs(sumET / GeV)))

            # transverse mass
            elET = Electron.fourvect.Pt()
            tauET = Tau.fourvect.Pt()
            elPhiVector = Vector2(Electron.fourvect.Px(), Electron.fourvect.Py())
            tauPhiVector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            dPhi_MET_electron = elPhiVector.DeltaPhi(MET_vect)
            dPhi_MET_tau  = tauPhiVector.DeltaPhi(MET_vect)
            mT = sqrt(2*MET*elET*(1 - cos(dPhi_MET_electron)))
            mTtau = sqrt(2*MET*tauET*(1 - cos(dPhi_MET_tau)))
            tree.mass_transverse_met_electron = mT
            tree.mass_transverse_met_tau = mTtau
            tree.dphi_met_electron = dPhi_MET_electron

            # ddR
            tree.ddr_tau_electron, tree.dr_tau_electron, tree.higgs_pt = eventshapes.DeltaDeltaR(Tau.fourvect, Electron.fourvect, MET_vect)
 
            """
            Higgs fancier mass calculation
            """
            collin_mass, tau_x, electron_x = mass.collinearmass(Tau, Electron, METx, METy)
            tree.mass_collinear_tau_electron = collin_mass
            tree.tau_x = tau_x
            tree.electron_x = electron_x
            mmc_mass, mmc_pt, mmc_met = mass.missingmass(Tau, Electron, METx, METy, sumET, 1)
            tree.mass_mmc_tau_electron = mmc_mass
            tree.pt_mmc_tau_electron = mmc_pt
            tree.met_mmc_tau_electron = mmc_met


            """
            Calculate fancier quantities for BDT input
            """

            tau2Vector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            electron2Vector = Vector2(Electron.fourvect.Px(), Electron.fourvect.Py()) # <--- is the electron cluster Et multiplied by (cos(phi),sin(phi)) where phi is obtained as recommended by EGamma group

            tree.met_phi_centrality = eventshapes.phi_centrality(tau2Vector, electron2Vector, MET_vect)

            mass_j1_j2 = -1111
            eta_product_j1_j2 = -1111
            eta_delta_j1_j2 = -1111
            tau_centrality_j1_j2 = -1111
            electron_centrality_j1_j2 = -1111
            tau_j1_j2_phi_centrality = -1111
            sphericity = -1111
            aplanarity = -1111

            #Calculate effective number of objects and mass of jets
            PtSum = 0
            PtSum2 = 0

            allJets = LorentzVector()

            for jet in event.jets:
                PtSum  += jet.fourvect.Pt()
                PtSum2 += (jet.fourvect.Pt())**2
                allJets += jet.fourvect

            leadJetPt = 0.0
            if len(event.jets) > 0:
                leadJetPt = event.jets[0].fourvect.Pt()
            tree.leadJetPt = leadJetPt

            PtSum  += (Tau.fourvect.Pt() + Electron.fourvect.Pt())
            PtSum2 += (Tau.fourvect.Pt()**2 + Electron.fourvect.Pt()**2)

            tree.neff_pt = PtSum2/(PtSum**2)
            tree.mass_all_jets = allJets.M()

            if len(event.jets) >= 2:
                event.jets.sort(key=lambda jet: jet.pt, reverse=True)

                jet1 = event.jets[0].fourvect
                jet2 = event.jets[1].fourvect

                jet1_2Vector = Vector2(jet1.Px(), jet1.Py())
                jet2_2Vector = Vector2(jet2.Px(), jet2.Py())

                tau_j1_j2_phi_centrality = eventshapes.phi_centrality(jet1_2Vector, jet2_2Vector, tau2Vector)

                mass_j1_j2 = (jet1 + jet2).M()
                eta_product_j1_j2 = jet1.Eta() * jet2.Eta()
                eta_delta_j1_j2 = abs(jet1.Eta() - jet2.Eta())
                tau_centrality_j1_j2 = eventshapes.eta_centrality(Tau.fourvect.Eta(), jet1.Eta(), jet2.Eta())
                electron_centrality_j1_j2 = eventshapes.eta_centrality(Electron.fourvect.Eta(), jet1.Eta(), jet2.Eta())

                sphericity, aplanarity = eventshapes.sphericity_aplanarity([Tau.fourvect,
                                                                            Electron.fourvect,
                                                                            jet1,
                                                                            jet2])


            tree.mass_j1_j2 = mass_j1_j2
            tree.eta_product_j1_j2 = eta_product_j1_j2
            tree.eta_delta_j1_j2 = eta_delta_j1_j2
            tree.tau_centrality_j1_j2 = tau_centrality_j1_j2
            tree.electron_centrality_j1_j2 = electron_centrality_j1_j2
            tree.tau_j1_j2_phi_centrality = tau_j1_j2_phi_centrality

            tree.sphericity = sphericity
            tree.aplanarity = aplanarity

            tree.nvtx = event.vxp_n

            event_weight = 1.0

            # Get pileup weight
            if self.metadata.datatype == datasets.MC:
                # set the event pileup weight
                event_weight = pileup_tool.GetCombinedWeight(event.RunNumber,
                                                             event.mc_channel_number,
                                                             event.averageIntPerXing)

                # Correct for trigger luminosity in period I-K
                if 185353 <= event.RunNumber <= 187815 and not event.EF_e22_medium:
                    event_weight *= 0.49528

                # Tau/Electron misidentification correction
                event_weight *= TauEfficiencySF(event, self.metadata.datatype)

                # Electron scale factors
                event_weight *= ElectronSF(event, self.metadata.datatype, pileup_tool)

                # ggF Reweighting
                event_weight *= ggFreweighting(event, self.metadata.name)

            tree.weight = event_weight

            # fill output ntuple
            tree.cutflow = cutflow.int()
            tree.Fill()

        self.output.cd()
        tree_train.FlushBaskets()
        tree_train.Write()
        tree_test.FlushBaskets()
        tree_test.Write()

        if self.metadata.datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        merged_cutflow.Write()

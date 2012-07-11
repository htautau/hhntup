import ROOT
import math

from argparse import ArgumentParser

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
from higgstautau import mass
from higgstautau.trigger import utils as triggerutils
from higgstautau.pileup import TPileupReweighting
from higgstautau.systematics import Systematics
from higgstautau.jetcalibration import JetCalibration

from goodruns import GRL
import subprocess

import random

YEAR = 2011

class muLHProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def __init__(self, options, **kwargs):

        super(muLHProcessor, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-type', default='None')
        parser.add_argument('--syst-term', default='None')
        self.args = parser.parse_args(options)
        self.args.syst_type = eval(self.args.syst_type)
        self.args.syst_term = eval(self.args.syst_term)
        

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
        #trigger_config = get_trigger_config()

        OutputModel = RecoTauLepBlock + EventVariables + RecoMET

        onfilechange = []
        if self.metadata.datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file, tree):
                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))

        # update the trigger config maps on every file change
        # onfilechange.append((update_trigger_config, (trigger_config,)))

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
        tree_train = Tree(name='lh_train', model=OutputModel)
        tree_test = Tree(name='lh_test', model=OutputModel)


        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']

        tree_train.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)
        tree_test.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)

        chain.always_read(copied_variables)

        # set the event filters

        Trigger = muMCSLTriggers
        if self.metadata.datatype == datasets.DATA:
            Trigger = muSLTriggers
        if self.metadata.datatype == datasets.EMBED:
            Trigger = noTriggers

        

        
        # passthrough for MC for trigger acceptance studies
        event_filters = EventFilterList([
            SetElectronsFourVector(),
            Trigger(),
            GRLFilter(self.grl, passthrough=self.metadata.datatype != datasets.DATA),
            JetCalibration(
                year=YEAR,
                datatype=self.metadata.datatype,
                verbose=False),
            PriVertex(),
            MuonPtSmearing(datatype=self.metadata.datatype),
            EgammaERescaling(datatype=self.metadata.datatype),
            Systematics(
                systematic_type=self.args.syst_type,
                systematic_term=self.args.syst_term,
                year=YEAR,
                datatype=self.metadata.datatype,
                verbose=False),
            JetPreSelection(),
            MuonPreSelection(),
            ElectronPreSelection(),
            JetOverlapRemoval(),
            JetCleaning(eta_max = 9999.0),
            ElectronLArHole(),
            TauLArHole(),
            LArHole(datatype=self.metadata.datatype),
            LArError(),
            LeptonOverlapRemoval(),
            DileptonVeto(),
            MuonSelection(),
            TauPreSelection(),
            TauSelection(),
            JetSelection(),
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
        chain.define_collection(name="truthjets", prefix="jet_antikt4truth_", size="jet_antikt4truth_n", mix=FourMomentum)
        chain.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        # define tree objects
        tree_train.define_object(name='tau', prefix='tau_')
        tree_train.define_object(name='lep', prefix='lep_')
        tree_test.define_object(name='tau', prefix='tau_')
        tree_test.define_object(name='lep', prefix='lep_')

        if self.metadata.datatype == datasets.MC:
            # Initialize the pileup reweighting tool
            pileup_tool = TPileupReweighting()
            pileup_tool.AddConfigFile('higgstautau/pileup/%s_defaults.prw.root' % self.metadata.category)
            pileup_tool.AddLumiCalcFile('grl/2011/lumicalc/lephad/ilumicalc_histograms_None_178044-191933.root')
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetUnrepresentedDataAction(2)
            pileup_tool.Initialize()
            print pileup_tool.getIntegratedLumiVector()


        # entering the main event loop...
        for event in chain:

            tree_train.reset()
            tree_test.reset()
            cutflow.reset()

            #Select if the event goes into the training or the testing tree
            tree = None
            if random.random() < 0.5:
                tree = tree_train
            else:
                tree = tree_test

            # Select tau with highest BDT score and surviving muon
            Tau = event.taus[0]
            Muon = event.muons[0]

            """
            RecoTauLepBlock filling
            """
            RecoTauLepBlock.set(event, tree, Tau, Muon, 0, ( self.metadata.datatype == datasets.MC ))


            """
            Jets
            """
            numJets = len(event.jets)
            tree.numJets = numJets

            numJets50 = 0
            numJets35 = 0

            for jet in event.jets:
                tree.jet_fourvect.push_back(jet.fourvect)
                tree.jet_jvtxf.push_back(jet.jvtxf)
                tree.jet_btag.push_back(jet.flavor_weight_JetFitterCOMBNN)
                if jet.fourvect.Pt() > 35*GeV:
                    numJets35 += 1
                    if jet.fourvect.Pt() > 50*GeV:
                        numJets50 += 1

            tree.numJets50 = numJets50
            tree.numJets35 = numJets35


            """
            Truth jets
            """

            if ( self.metadata.datatype == datasets.MC ):
                for truthjet in event.truthjets:
                    tree.truthjet_fourvect.push_back(truthjet.fourvect)
                


            """
            Miscellaneous
            """
            tree.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or
                                         (vtx.type == 3 and vtx.nTracks >= 2)])

            sumPt = 0

            sumPt += Tau.fourvect.Pt()
            sumPt += Muon.fourvect.Pt()
            for jet in event.jets:
                sumPt += jet.fourvect.Pt()

            tree.sumPt = sumPt

            #missing ET
            METx = event.MET_RefFinal_BDTMedium_etx
            METy = event.MET_RefFinal_BDTMedium_ety
            MET_vect = Vector2(METx, METy)
            MET = MET_vect.Mod()
            tree.MET = MET
            getattr(tree, 'MET_vect').set_from(MET_vect)
            sumET = event.MET_RefFinal_BDTMedium_sumet
            tree.MET_sig = (2. * MET_vect.Mod() / GeV) / (utils.sign(sumET) * sqrt(abs(sumET / GeV)))

            #transverse mass
            muET = Muon.fourvect.Pt()
            tauET = Tau.fourvect.Pt()
            muPhiVector = Vector2(Muon.fourvect.Px(), Muon.fourvect.Py())
            tauPhiVector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            dPhi_MET_muon = muPhiVector.DeltaPhi(MET_vect)
            dPhi_MET_tau  = tauPhiVector.DeltaPhi(MET_vect)
            mT = sqrt(2*MET*muET*(1 - cos(dPhi_MET_muon)))
            mTtau = sqrt(2*MET*tauET*(1 - cos(dPhi_MET_tau)))
            tree.mass_transverse_met_lep = mT
            tree.mass_transverse_met_tau = mTtau
            tree.dphi_met_lep = dPhi_MET_muon

            #ddR
            tree.ddr_tau_lep, tree.dr_tau_lep, tree.resonance_pt_tau_lep = eventshapes.DeltaDeltaR(Tau.fourvect, Muon.fourvect, MET_vect)
            

            """
            Higgs fancier mass calculation
            """
            collin_mass, tau_x, muon_x = mass.collinearmass(Tau, Muon, METx, METy)
            tree.mass_collinear_tau_lep = collin_mass
            tree.tau_x  = tau_x
            tree.lep_x = muon_x
            mmc_mass, mmc_pt, mmc_met = mass.missingmass(Tau, Muon, METx, METy, sumET, 1)
            tree.mass_mmc_tau_lep = mmc_mass
            tree.pt_mmc_tau_lep = mmc_pt
            tree.met_mmc_tau_lep = mmc_met



            """
            Calculate fancier quantities for BDT input
            """

            tau2Vector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            muon2Vector = Vector2(Muon.fourvect.Px(), Muon.fourvect.Py())

            tree.met_phi_centrality = eventshapes.phi_centrality(tau2Vector, muon2Vector, MET_vect)

            mass_j1_j2 = -1111
            eta_product_j1_j2 = -1111
            eta_delta_j1_j2 = -1111
            tau_centrality_j1_j2 = -1111
            muon_centrality_j1_j2 = -1111
            tau_j1_j2_phi_centrality = -1111
            sphericity = -1111
            aplanarity = -1111

            #Calculate jet stuff
            allJetList = []

            for jet in event.jets:
                allJetList.append(jet.fourvect)

            leadJetPt = 0.0
            if len(event.jets) > 0:
                leadJetPt = event.jets[0].fourvect.Pt()
            tree.leadJetPt = leadJetPt

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
                muon_centrality_j1_j2 = eventshapes.eta_centrality(Muon.fourvect.Eta(), jet1.Eta(), jet2.Eta())

            sphericity, aplanarity = eventshapes.sphericity_aplanarity([Tau.fourvect, Muon.fourvect] + allJetList)


            tree.mass_j1_j2 = mass_j1_j2
            tree.eta_product_j1_j2 = eta_product_j1_j2
            tree.eta_delta_j1_j2 = eta_delta_j1_j2
            tree.tau_centrality_j1_j2 = tau_centrality_j1_j2
            tree.lep_centrality_j1_j2 = muon_centrality_j1_j2
            tree.tau_j1_j2_phi_centrality = tau_j1_j2_phi_centrality

            tree.sphericity = sphericity
            tree.aplanarity = aplanarity

            tree.nvtx = event.vxp_n

            event_weight = 1.0

            #Get pileup weight
            if self.metadata.datatype == datasets.MC:
                # set the event pileup weight
                event_weight = pileup_tool.GetCombinedWeight(event.RunNumber,
                                                             event.mc_channel_number,
                                                             event.averageIntPerXing)

                # Correct for trigger luminosity in period I-K
                if 185353 <= event.RunNumber <= 187815 and not event.EF_mu18_MG_medium:
                    event_weight *= 0.29186

                #Tau/Electron misidentification correction
                event_weight *= TauEfficiencySF(event, self.metadata.datatype)

                #Muon Scale factors
                event_weight *= MuonSF(event, self.metadata.datatype, pileup_tool)

                #ggF Reweighting
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

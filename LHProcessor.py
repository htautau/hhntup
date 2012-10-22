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
from higgstautau.corrections import reweight_ggf
from higgstautau.embedding import EmbeddingPileupPatch

from goodruns import GRL
import subprocess

YEAR = 2011
VERBOSE = False

class LHProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def __init__(self, options, **kwargs):

        super(LHProcessor, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-terms', default=None)
        self.args = parser.parse_args(options)
        if self.args.syst_terms is not None:
            self.args.syst_terms = [
                eval('Systematics.%s' % term) for term in
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

        if self.metadata.datatype != datasets.EMBED:
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

        ## set the trigger filter
        Trigger = noTriggers
        if self.metadata.datatype == datasets.DATA:
            if self.metadata.stream == 'Egamma':
                Trigger = eSLTriggers
            if self.metadata.stream == 'Muons':
                Trigger = muSLTriggers
            if self.metadata.stream == 'JetTauEtmiss':
                Trigger = AllDataLTTriggers

        if self.metadata.datatype == datasets.MC:
            Trigger = AllMCTriggers


        ## Setting event filters
        event_filters = EventFilterList([
            PrepareInputTree(),
            Trigger(),
            GRLFilter( self.grl, passthrough=self.metadata.datatype == datasets.MC ),
            #EmbeddingPileupPatch( passthrough=self.metadata.datatype != datasets.EMBED ),
            JetCalibration( year=YEAR, datatype=self.metadata.datatype, verbose=False ),
            PriVertex(),
            Systematics( terms=self.args.syst_terms, year=YEAR, datatype=self.metadata.datatype, verbose=VERBOSE ),
            JetPreSelection(),
            MuonPtSmearing( datatype=self.metadata.datatype ),
            MuonPreSelection(),
            EgammaERescaling( datatype=self.metadata.datatype ),
            ElectronPreSelection(),
            JetOverlapRemoval(),
            ElectronEtaSelection(),
            JetCleaning( self.metadata.datatype, YEAR ),
            ElectronLArHole(),
            TauHasTrack(1),
            TauLArHole(1),
            LArHole( datatype=self.metadata.datatype ),
            LArError(),
            LeptonOverlapRemoval(),
            DileptonVeto(),
            LeptonSelection( datatype=self.metadata.datatype, stream=self.metadata.stream ),
            TauPreSelection(),
            TauSelection(),
            JetSelection(),
            FinalOverlapRemoval(),
            ElectronIsoCorrection( datatype=self.metadata.datatype ),
        #AntiVBFFilter()
        ])

        self.filters['event'] = event_filters

        chain.filters += event_filters

        if self.metadata.datatype != datasets.EMBED:
            cutflow = Cutflow()

        # define tree collections
        chain.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n", mix=FourMomentum)
        chain.define_collection(name="electrons", prefix="el_", size="el_n", mix=ElFourMomentum)
        chain.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        chain.define_collection(name="truthjets", prefix="jet_antikt4truth_", size="jet_antikt4truth_n", mix=FourMomentum)
        chain.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")
        chain.define_object(name='isLTT', prefix='')
        chain.define_object(name='leptonType', prefix='')

        # define tree objects
        tree_train.define_object(name='tau', prefix='tau_')
        tree_train.define_object(name='lep', prefix='lep_')
        tree_test.define_object(name='tau', prefix='tau_')
        tree_test.define_object(name='lep', prefix='lep_')

        if self.metadata.datatype == datasets.MC or self.metadata.datatype == datasets.EMBED:
            from externaltools import PileupReweighting
            from ROOT import Root
            # Initialize the pileup reweighting tool
            pileup_tool = Root.TPileupReweighting()
            if YEAR == 2011:
                pileup_tool.AddConfigFile(PileupReweighting.get_resource('mc11b_defaults.prw.root'))
                pileup_tool.AddLumiCalcFile('lumi/2011/hadhad/ilumicalc_histograms_None_178044-191933.root')
            elif YEAR == 2012:
                pileup_tool.AddConfigFile(PileupReweighting.get_resource('mc12a_defaults.prw.root'))
                pileup_tool.SetDataScaleFactors(1./1.11)
                pileup_tool.AddLumiCalcFile('lumi/2012/hadhad/ilumicalc_histograms_None_200841-205113.root')
            else:
                raise ValueError('No pileup reweighting defined for year %d' %
                        YEAR)
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetUnrepresentedDataAction(2)
            pileup_tool.Initialize()
            print pileup_tool.getIntegratedLumiVector()


        # entering the main event loop...
        for event in chain:

            #tree_train.reset() Use reset=True in the Fill below
            #tree_test.reset()
            if self.metadata.datatype != datasets.EMBED:
                cutflow.reset()

            #Select if the event goes into the training or the testing tree
            tree = None
            if event.EventNumber % 2:
                tree = tree_train
            else:
                tree = tree_test

            # Select tau with highest BDT score and surviving lepton
            try:
                Tau = event.taus[0]
                Lep = None
                leptype = None
                if event.leptonType == 'mu':
                    Lep = event.muons[0]
                    leptype = 0
                if event.leptonType == 'e':
                    Lep = event.electrons[0]
                    leptype = 1
            except IndexError:
                continue

            """
            RecoTauLepBlock filling
            """
            RecoTauLepBlock.set(event, tree, Tau, Lep, leptype, ( self.metadata.datatype == datasets.MC ))


            """
            Jets
            """
            numJets = len(event.jets)
            tree.numJets = numJets

            numJets50 = 0
            numJets30 = 0

            for jet in event.jets:
                tree.jet_fourvect.push_back(jet.fourvect)
                tree.jet_jvtxf.push_back(jet.jvtxf)
                tree.jet_btag.push_back(jet.flavor_weight_JetFitterCOMBNN)
                if jet.fourvect.Pt() > 30*GeV:
                    numJets30 += 1
                    if jet.fourvect.Pt() > 50*GeV:
                        numJets50 += 1

            tree.numJets50 = numJets50
            tree.numJets30 = numJets30


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
            sumPt += Lep.fourvect.Pt()
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
            lepET = Lep.fourvect.Pt()
            tauET = Tau.fourvect.Pt()
            lepPhiVector = Vector2(Lep.fourvect.Px(), Lep.fourvect.Py())
            tauPhiVector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            dPhi_MET_lep = lepPhiVector.DeltaPhi(MET_vect)
            dPhi_MET_tau  = tauPhiVector.DeltaPhi(MET_vect)
            mT = sqrt(2*MET*lepET*(1 - cos(dPhi_MET_lep)))
            mTtau = sqrt(2*MET*tauET*(1 - cos(dPhi_MET_tau)))
            tree.mass_transverse_met_lep = mT
            tree.mass_transverse_met_tau = mTtau
            tree.dphi_met_lep = dPhi_MET_lep

            #ddR
            tree.ddr_tau_lep, tree.dr_tau_lep, tree.resonance_pt_tau_lep = eventshapes.DeltaDeltaR(Tau.fourvect, Lep.fourvect, MET_vect)


            """
            Higgs fancier mass calculation
            """
            collin_mass, tau_x, lep_x = mass.collinearmass(Tau, Lep, METx, METy)
            tree.mass_collinear_tau_lep = collin_mass
            tree.tau_x  = tau_x
            tree.lep_x = lep_x
            mmc_mass, mmc_pt, mmc_met = mass.missingmass(Tau, Lep, METx, METy, sumET, leptype)
            tree.mass_mmc_tau_lep = mmc_mass
            tree.pt_mmc_tau_lep = mmc_pt.Pt()
            tree.met_mmc_tau_lep = mmc_met.Mod()



            """
            Calculate fancier quantities for BDT input
            """

            tau2Vector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            lep2Vector = Vector2(Lep.fourvect.Px(), Lep.fourvect.Py())

            tree.met_phi_centrality = eventshapes.phi_centrality(tau2Vector, lep2Vector, MET_vect)

            mass_j1_j2 = -1111
            eta_product_j1_j2 = -1111
            eta_delta_j1_j2 = -1111
            tau_centrality_j1_j2 = -1111
            lep_centrality_j1_j2 = -1111
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
                lep_centrality_j1_j2 = eventshapes.eta_centrality(Lep.fourvect.Eta(), jet1.Eta(), jet2.Eta())

            sphericity, aplanarity = eventshapes.sphericity_aplanarity([Tau.fourvect, Lep.fourvect] + allJetList)


            tree.mass_j1_j2 = mass_j1_j2
            tree.eta_product_j1_j2 = eta_product_j1_j2
            tree.eta_delta_j1_j2 = eta_delta_j1_j2
            tree.tau_centrality_j1_j2 = tau_centrality_j1_j2
            tree.lep_centrality_j1_j2 = lep_centrality_j1_j2
            tree.tau_j1_j2_phi_centrality = tau_j1_j2_phi_centrality

            tree.sphericity = sphericity
            tree.aplanarity = aplanarity

            tree.nvtx = event.vxp_n

            event_weight = 1.0

            npileup_vtx = len([vtx for vtx in event.vertices
                           if pileup_vertex_selection(vtx)])

            ## Event weight corrections for MC
            ###################################
            if self.metadata.datatype == datasets.MC:
                ## Apply mc event weight:
                try:
                    event_weight *= event.mcevt_weight[0][0]
                except AttributeError:
                    pass

                # set the event pileup weight
                event_weight = pileup_tool.GetCombinedWeight(event.RunNumber,
                                                             event.mc_channel_number,
                                                             event.averageIntPerXing)

                # Correct for trigger luminosity in period I-K
                if event.isLTT:
                    if 185353 <= event.RunNumber <= 187815 and not event.EF_mu18_MG_medium:
                        event_weight *= 0.49528
                else:
                    if 185353 <= event.RunNumber <= 187815 and not event.EF_mu18_MG_medium:
                        event_weight *= 0.29186

                #Tau/Electron misidentification correction
                event_weight *= TauEfficiencySF(event, self.metadata.datatype)

                #Lepton Efficiency scale factors
                if event.leptonType == 'mu':
                    event_weight *= MuonSF(event, self.metadata.datatype, pileup_tool)
                if event.leptonType == 'e':
                    event_weight *= ElectronSF(event, self.metadata.datatype)

                #Lepton Trigger scale factors
                if not event.isLTT:
                    event_weight *= LeptonSLTSF(event, self.metadata.datatype)
                else:
                    if event.leptonType == 'mu':
                        event_weight *= MuonLTTSF(Lep, event.RunNumber, self.metadata.datatype)
                    if event.leptonType == 'e':
                        event_weight *= ElectronLTTSF(Lep, event.RunNumber, self.metadata.datatype)

                #ggF Reweighting
                event_weight *= reweight_ggf(event, self.metadata.name)


             ## Event weight corrections for embedded samples
            #################################################
            if self.metadata.datatype == datasets.EMBED:

                mc_w              = 1.0
                leptonsf_w        = 1.0
                leptontriggersf_w = 1.0
                tautriggersf_w    = 1.0

                ## Apply mc event weight:
                try:
                    mc_w = event.mcevt_weight[0][0]
                except AttributeError:
                    pass

                #Lepton Efficiency scale factors
                if event.leptonType == 'mu':
                    leptonsf_w = MuonSF(event, self.metadata.datatype, pileup_tool)
                if event.leptonType == 'e':
                    leptonsf_w = ElectronSF(event, self.metadata.datatype)

                #Lepton Trigger scale factors
                if not event.isLTT:
                    leptontriggersf_w = LeptonSLTSF(event, self.metadata.datatype)
                else:
                    if event.leptonType == 'mu':
                        leptontriggersf_w = MuonLTTSF(Lep, event.RunNumber, self.metadata.datatype)
                    if event.leptonType == 'e':
                        leptontriggersf_w = ElectronLTTSF(Lep, event.RunNumber, self.metadata.datatype)

                    #Tau trigger scale factors
                    tautriggersf_w = EmbedTauTriggerCorr(Tau, npileup_vtx, event.RunNumber)

                # print '---------------------'
                # print 'event is LTT : %s' % str(event.isLTT)
                # print 'leptype      : %s' % event.leptonType
                # print 'tauTriggerSF : %f' % tautriggersf_w
                # print 'leptonSF     : %f' % leptonsf_w
                # print 'lepTriggerSF : %f' % leptontriggersf_w

                event_weight = mc_w * tautriggersf_w * leptonsf_w * leptontriggersf_w

            tree.weight = event_weight

            # fill output ntuple
            if self.metadata.datatype != datasets.EMBED:
                tree.cutflow = cutflow.int()
            tree.Fill(reset=True)

        self.output.cd()
        tree_train.FlushBaskets()
        tree_train.Write()
        tree_test.FlushBaskets()
        tree_test.Write()

        if self.metadata.datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        if self.metadata.datatype != datasets.EMBED:
            merged_cutflow.Write()

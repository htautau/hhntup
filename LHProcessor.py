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

        ## Use the metadata to set the year
        YEAR = self.metadata.year
        
        ## Initiate the MMC object
        mmc = mass.MMC(year=YEAR, channel='lh')

        ## Set the output ntuple model
        OutputModel = RecoTauLepBlock + EventVariables + RecoMET + SysWeights

        onfilechange = []
        if self.metadata.datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file, tree):
                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))


        if self.metadata.datatype != datasets.EMBED:
            merged_cutflow = Hist(7, 0, 7, name='cutflow', type='D')

            def update_cutflow(student, cutflow, name, file, tree):
                if year == 2011:
                    cutflow += file.cutflow
                elif datatype == datasets.MC:
                    cutflow[0] += file.cutflow_event[0]
                    cutflow[1] += file.cutflow_event_mc_weight[0]
                else:
                    cutflow[0] += file.cutflow_event[0]

                onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        chain = TreeChain(self.metadata.treename,
                         files=self.files,
                         events=self.events,
                         cache=True,
                         cache_size=10000000,
                         learn_entries=30,
                         onfilechange=onfilechange,
                         verbose=True)

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


        #####################################################################################
        ## Pileup tool for pileup reweighting (and so much more)
        
        if self.metadata.datatype == datasets.MC or self.metadata.datatype == datasets.EMBED:
            from externaltools import PileupReweighting
            from ROOT import Root
            # Initialize the pileup reweighting tool
            pileup_tool = Root.TPileupReweighting('AdvancedPURT')
            if YEAR == 2011:
                pileup_tool.AddConfigFile(PileupReweighting.get_resource('mc11b_defaults.root'))
                pileup_tool.SetDataScaleFactors(1.0)
                pileup_tool.AddLumiCalcFile('lumi/2011/hadhad/ilumicalc_histograms_None_178044-191933.root')
            elif YEAR == 2012:
                pileup_tool.AddConfigFile(PileupReweighting.get_resource('mc12a_defaults.prw.root'))
                pileup_tool.SetDataScaleFactors(1./1.11)
                pileup_tool.AddLumiCalcFile('lumi/2012/hadhad/ilumicalc_histograms_None_200841-205113.root')
            else:
                raise ValueError('No pileup reweighting defined for year %d' %
                        YEAR)
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetDefaultChannel(0)
            pileup_tool.SetUnrepresentedDataAction(2)
            pileup_tool.Initialize()

        
        ######################################################################################
        ## Instantiate correction tools
        
        ## Muon momentum corrections
        MMCorrections = None
        from externaltools import MuonMomentumCorrections as MMC
        from ROOT import MuonSmear

        if YEAR == 2011:
            MMCorrections = MuonSmear.SmearingClass('Data11',
                                                    'staco',
                                                    'pT',
                                                    'Rel17',
                                                    MMC.RESOURCE_PATH)
            MMCorrections.UseScale(1)
            MMCorrections.UseImprovedCombine()
            MMCorrections.RestrictCurvatureCorrections(2.5)
            MMCorrections.FillScales('KC')

        if YEAR == 2012:
            MMCorrections = MuonSmear.SmearingClass('Data12',
                                                    'staco',
                                                    'pT',
                                                    'Rel17.2_preliminary',
                                                    MMC.RESOURCE_PATH)
            MMCorrections.UseScale(1)
            MMCorrections.UseImprovedCombine()


        ## Muon isolation corrections
        MICorrections = None
        from externaltools.bundle_2012 import MuonIsolationCorrection
        from ROOT import CorrectCaloIso
        MICorrections = CorrectCaloIso()


        ## Egamma momentum corrections
        EGMCorrections = None
        if YEAR == 2011:
            from externaltools import egammaAnalysisUtils
            from ROOT import eg2011
            EGMCorrections = eg2011.EnergyRescaler()

        if YEAR == 2012:
            from externaltools import egammaAnalysisUtils
            from ROOT import eg2011
            EGMCorrections = eg2011.EnergyRescaler()


        ## Egamma isolation corrections
        from ROOT import CaloIsoCorrection
        EGICorrections = CaloIsoCorrection


        ## Muon isolation efficiency correction
        from rootpy.io import open as ropen
        HERE = os.path.dirname(os.path.abspath(__file__))
        MuonIsoCorrFile = ropen(os.path.join(HERE, 'higgstautau/lephad/SF_2D_UptoE5.root'))
        MIEffCorrections = MuonIsoCorrFile.Get('SF_2D_PtVsNvx')


        if self.metadata.datatype != datasets.DATA:
            ## Muon Efficiency corrections
            MEffCorrections = None
            from externaltools import MuonEfficiencyCorrections
            from ROOT import Analysis

            if YEAR == 2011:
                int_lum = pileup_tool.getIntegratedLumiVector()
                MEffCorrections = Analysis.AnalysisMuonEfficiencyScaleFactors('STACO_CB',
                                                                              int_lum,
                                                                              'MeV',
                                                                              MuonEfficiencyCorrections.RESOURCE_PATH)

            if YEAR == 2012:
                MEffCorrections = Analysis.AnalysisMuonConfigurableScaleFactors('',
                                                                                MuonEfficiencyCorrections.RESOURCE_PATH + '/STACO_CB_2012_SF.txt',
                                                                                'MeV',
                                                                                Analysis.AnalysisMuonConfigurableScaleFactors.AverageOverPeriods)
                MEffCorrections.setRunInterval(200804, 210308)
                MEffCorrections.Initialise()


            ## Electron Efficiency corrections
            EGEffCorrections = None
            from ROOT import egammaSFclass
            EGEffCorrections = egammaSFclass()

        
            ## LTT trigger corrections
            LTTCorrections = None
            if YEAR == 2011:
                from externaltools.bundle_2011 import HSG4TriggerSF as HSG4
            if YEAR == 2012:
                from externaltools.bundle_2012 import HSG4TriggerSF as HSG4
            from ROOT import HSG4TriggerSF
            LTTCorrections = HSG4TriggerSF(HSG4.RESOURCE_PATH)

            LTTtauCorrections = None
            if YEAR == 2011:
                from externaltools.bundle_2011 import TauTriggerCorrections  as TTC
                LTTtauPath = TTC.RESOURCE_PATH
            if YEAR == 2012:
                from externaltools.bundle_2012 import TauTriggerCorrections  as TTC
                LTTtauPath = TTC.RESOURCE_PATH

            from ROOT import TauTriggerCorrections
            LTTtauCorrections1 = TauTriggerCorrections()
            LTTtauCorrections2 = TauTriggerCorrections()
            LTTtauCorrections3 = TauTriggerCorrections()

            if YEAR == 2011:
                LTTtauCorrections1.loadInputFile(os.path.join(LTTtauPath, 'triggerSF_wmcpara_EF_tau20_medium1.root'), '1P3P', 'BDTm')
                LTTtauCorrections2.loadInputFile(os.path.join(LTTtauPath, 'triggerSF_EF_tau20_medium1.root'))
                LTTtauCorrections3.loadInputFile(os.path.join(LTTtauPath, 'triggerSF_EF_tau16_loose.root'))
            if YEAR == 2012:
                LTTtauCorrections1.loadInputFile(os.path.join(LTTtauPath, 'triggerSF_EF_tau20_medium1.root'))
                LTTtauCorrections2.loadInputFile(os.path.join(LTTtauPath, 'triggerSF_EF_tau20Ti_medium1.root'))


            ## SLT trigger corrections
            SLTCorrections = None
            if YEAR == 2011:
                from externaltools import TrigMuonEfficiency
                from ROOT import LeptonTriggerSF
                from ROOT import TrigMuonEff
                SLTCorrections = LeptonTriggerSF(TrigMuonEfficiency.RESOURCE_PATH)
            if YEAR == 2012:
                from externaltools import TrigMuonEfficiency
                from ROOT import LeptonTriggerSF
                from ROOT import TrigMuonEff
                SLTCorrections = LeptonTriggerSF(2012, TrigMuonEfficiency.RESOURCE_PATH, 'muon_trigger_sf_2012_AtoE.root')


        #####################################################################################
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
        #AcceptanceChallenge(),
            PrepareInputTree(),
            Trigger( year=YEAR ),
            pseudoGRLFilter(),
            GRLFilter( self.grl, passthrough=(self.metadata.datatype != datasets.DATA) ),
            #EmbeddingPileupPatch( passthrough=self.metadata.datatype != datasets.EMBED ),
            JetCalibration( year=YEAR, datatype=self.metadata.datatype, verbose=False ),
            MuonPtSmearing( datatype=self.metadata.datatype, year=YEAR, tool=MMCorrections ),
            EgammaERescaling( datatype=self.metadata.datatype, year=YEAR, tool=EGMCorrections ),
            Systematics( terms=self.args.syst_terms, year=YEAR, channel='lh', datatype=self.metadata.datatype, verbose=VERBOSE ),
            #METCorrection( datatype=self.metadata.datatype, year=YEAR ),
            PriVertex(),
            JetPreSelection(),
            MuonOverlapSelection(year=YEAR),
            TauMuonOverlapRemoval(),
            MuonPreSelection(year=YEAR),
            ElectronPreSelection(year=YEAR),
            JetOverlapRemoval(),
            ElectronEtaSelection(),
            JetCleaning( self.metadata.datatype, year=YEAR ),
            ElectronLArHole(),
            TauHasTrack(1),
            TauLArHole(1),
            LArHole( datatype=self.metadata.datatype ),
            LArError(),
            LeptonOverlapRemoval(),
            DileptonVeto(),
            LeptonSelection( datatype=self.metadata.datatype, stream=self.metadata.stream, year=YEAR ),
            TauPreSelection(),
            TauSelection(),
            JetSelection( year=YEAR, bunny_ear_protection=False ),
            FinalOverlapRemoval(),
            ElectronIsoCorrection( datatype=self.metadata.datatype, year=YEAR, tool=EGICorrections ),
            MuonIsoCorrection(datatype=self.metadata.datatype, year=YEAR, tool=MICorrections )
            #AntiVBFFilter()
        ])

        self.filters['event'] = event_filters

        chain.filters += event_filters

        if self.metadata.datatype != datasets.EMBED:
            cutflow = Cutflow()

        # define tree collections
        chain.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n", mix=FourMomentum)
        chain.define_collection(name="electrons", prefix="el_", size="el_n", mix=FourMomentum)
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
                mmc_leptype = None
                if event.leptonType == 'mu':
                    Lep = event.muons[0]
                    leptype = 0
                    mmc_leptype = 1
                if event.leptonType == 'e':
                    Lep = event.electrons[0]
                    leptype = 1
                    mmc_leptype = 0
            except IndexError:
                continue

            """
            RecoTauLepBlock filling
            """
            RecoTauLepBlock.set(event, tree, Tau, Lep, leptype, ( self.metadata.datatype == datasets.MC ), year=YEAR)


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
                if YEAR == 2011:
                    tree.jet_btag.push_back(jet.flavor_weight_JetFitterCOMBNN)
                if YEAR == 2012:
                    tree.jet_btag.push_back(jet.flavor_weight_MV1)
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
            if YEAR == 2011:
                METx = event.MET_RefFinal_BDTMedium_etx
                METy = event.MET_RefFinal_BDTMedium_ety
                sumET = event.MET_RefFinal_BDTMedium_sumet
            if YEAR == 2012:
                METx = event.MET_RefFinal_STVF_etx
                METy = event.MET_RefFinal_STVF_ety
                sumET = event.MET_RefFinal_STVF_sumet
            MET_vect = Vector2(METx, METy)
            MET = MET_vect.Mod()
            tree.MET = MET
            tree.sumET = sumET
            getattr(tree, 'MET_vect').set_from(MET_vect)
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
            
            mmc_mass, mmc_pt, mmc_met = -1,-1,-1
            
            mmc_mass, mmc_pt, mmc_met = mmc.mass(Tau, Lep, METx, METy, sumET, mmc_leptype, njets25=numJets)
            mmc_pt = mmc_pt.Pt()
            mmc_met = mmc_met.Mod()
            
            tree.mass_mmc_tau_lep = mmc_mass
            tree.pt_mmc_tau_lep = mmc_pt
            tree.met_mmc_tau_lep = mmc_met



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

            
            tree.w_mc           = 1.0
            tree.w_pileup       = 1.0
            tree.w_lumi         = 1.0
            tree.w_tauesf       = 1.0
            tree.w_tauidsf      = 1.0
            tree.w_leptonsf     = 1.0
            tree.w_leptontrigsf = 1.0
            tree.w_tautriggersf = 1.0
            tree.w_muonisosf    = 1.0
            tree.w_ggf          = 1.0
            
            ## Event weight corrections for MC
            ###################################
            if self.metadata.datatype == datasets.MC:

                mc_w           = 1.0
                pileup_w       = 1.0
                lumi_w         = 1.0
                tauesf_w       = 1.0
                tauidsf_w      = 1.0
                leptonsf_w     = 1.0
                leptontrigsf_w = 1.0
                tautriggersf_w = 1.0
                muonisosf_w    = 1.0
                ggf_w          = 1.0
                
                ## Apply mc event weight:
                try:
                    mc_w = event.mcevt_weight[0][0]
                except AttributeError:
                    pass

                # set the event pileup weight

                pileup_w = pileup_tool.GetCombinedWeight(event.RunNumber,
                                                         0,
                                                         event.averageIntPerXing)

                # Correct for trigger luminosity in period I-K
                if YEAR == 2011:
                    if event.leptonType == 'mu':
                        if event.isLTT:
                            if 185353 <= event.RunNumber <= 187815 and not event.EF_tau20_medium_mu15:
                                lumi_w = 0.49528
                        else:
                            if 185353 <= event.RunNumber <= 187815 and not event.EF_mu18_MG_medium:
                                lumi_w = 0.29186
                    if event.leptonType == 'e':
                        if event.isLTT:
                            if 185353 <= event.RunNumber <= 187815 and not event.EF_tau20_medium_e15_medium:
                                lumi_w = 0.49528
                        else:
                            if 185353 <= event.RunNumber <= 187815 and not event.EF_e22_medium:
                                lumi_w = 0.49528

                # #Tau/Electron misidentification correction
                # tauesf_w, tree.sys_tau_ESF_UP, tree.sys_tau_ESF_DOWN = TauEfficiencySF(event,
                #                                                                        self.metadata.datatype,
                #                                                                        YEAR)

                #Tau ID scale factor correction
                # tauidsf_w, tree.sys_tau_IDSF_UP, tree.sys_tau_IDSF_DOWN = TauIDSF(event,
                #                                                                   self.metadata.datatype,
                #                                                                   YEAR)

                #Lepton Efficiency scale factors
                if event.leptonType == 'mu':
                    leptonsf_w, tree.sys_mu_EFFSF_UP, tree.sys_mu_EFFSF_DOWN = MuonSF(MEffCorrections,
                                                                                      SLTCorrections,
                                                                                      event,
                                                                                      self.metadata.datatype,
                                                                                      pileup_tool,
                                                                                      YEAR, event.RunNumber, event.isLTT, TrigMuonEff)
                if event.leptonType == 'e':
                    leptonsf_w, tree.sys_e_EFFSF_UP, tree.sys_e_EFFSF_DOWN = ElectronSF(EGEffCorrections,
                                                                                        event,
                                                                                        self.metadata.datatype,
                                                                                        pileup_tool,
                                                                                        YEAR, event.RunNumber, event.isLTT)

                #Lepton Trigger scale factors
                if not event.isLTT:
                    if event.leptonType == 'mu':
                        leptontrigsf_w, tree.sys_mu_TRIGSF_UP, tree.sys_mu_TRIGSF_DOWN = LeptonSLTSF(SLTCorrections,
                                                                                                     event,
                                                                                                     self.metadata.datatype,
                                                                                                     pileup_tool, 
                                                                                                     YEAR)
                else:
                    if event.leptonType == 'mu':
                        leptontrigsf_w, tree.sys_mu_TRIGSF_UP, tree.sys_mu_TRIGSF_DOWN = MuonLTTSF(LTTCorrections,
                                                                                                   event,
                                                                                                   self.metadata.datatype,
                                                                                                   YEAR,
                                                                                                   pileup_tool)
                    if event.leptonType == 'e':
                        leptontrigsf_w, tree.sys_e_TRIGSF_UP, tree.sys_e_TRIGSF_DOWN = ElectronLTTSF(LTTCorrections,
                                                                                                     event,
                                                                                                     self.metadata.datatype,
                                                                                                     YEAR,
                                                                                                     pileup_tool,
                                                                                                     event.RunNumber)
                        
                    #Tau trigger scale factors
                    tautriggersf_w, tree.sys_tau_TRIGSF_UP, tree.sys_tau_TRIGSF_DOWN = TauTriggerSF(LTTtauCorrections1,
                                                                                                    LTTtauCorrections2,
                                                                                                    LTTtauCorrections3,
                                                                                                    Tau,
                                                                                                    npileup_vtx,
                                                                                                    event.RunNumber,
                                                                                                    YEAR,
                                                                                                    event.leptonType,
                                                                                                    pileup_tool,
                                                                                                    self.metadata.datatype)

                if event.leptonType == 'mu':
                    muonisosf_w, tree.sys_mu_ISOSF_UP, tree.sys_mu_ISOSF_DOWN = MuonIsoEffCorrection(MIEffCorrections,
                                                                                                     event,
                                                                                                     self.metadata.datatype,
                                                                                                     YEAR)

                #ggF Reweighting
                ggf_w = reweight_ggf(event,
                                     self.metadata.name)

                event_weight = mc_w * pileup_w * lumi_w * tauesf_w * tauidsf_w * leptonsf_w * leptontrigsf_w * tautriggersf_w * muonisosf_w * ggf_w

                tree.w_mc           = mc_w
                tree.w_pileup       = pileup_w
                tree.w_lumi         = lumi_w
                tree.w_tauesf       = tauesf_w
                tree.w_tauidsf      = tauidsf_w
                tree.w_leptonsf     = leptonsf_w
                tree.w_leptontrigsf = leptontrigsf_w
                tree.w_tautriggersf = tautriggersf_w
                tree.w_muonisosf    = muonisosf_w
                tree.w_ggf          = ggf_w


             ## Event weight corrections for embedded samples
            #################################################
            if self.metadata.datatype == datasets.EMBED:

                mc_w              = 1.0
                leptonsf_w        = 1.0
                leptontrigsf_w    = 1.0
                tautriggersf_w    = 1.0
                tauidsf_w         = 1.0
                muonisosf_w       = 1.0

                ## Apply mc event weight:
                try:
                    mc_w = event.mcevt_weight[0][0]
                except AttributeError:
                    pass

                # #Tau/Electron misidentification correction
                # tauesf_w, tree.sys_tau_ESF_UP, tree.sys_tau_ESF_DOWN = TauEfficiencySF(event,
                #                                                                        self.metadata.datatype,
                #                                                                        YEAR)

                # #Tau ID scale factor correction
                # tauidsf_w, tree.sys_tau_IDSF_UP, tree.sys_tau_IDSF_DOWN = TauIDSF(event,
                #                                                                   self.metadata.datatype,
                #                                                                   YEAR)

                #Lepton Efficiency scale factors
                if event.leptonType == 'mu':
                    leptonsf_w, tree.sys_mu_EFFSF_UP, tree.sys_mu_EFFSF_DOWN = MuonSF(MEffCorrections,
                                                                                      SLTCorrections,
                                                                                      event,
                                                                                      self.metadata.datatype,
                                                                                      pileup_tool,
                                                                                      YEAR, event.RunNumber, event.isLTT, TrigMuonEff)
                if event.leptonType == 'e':
                    leptonsf_w, tree.sys_e_EFFSF_UP, tree.sys_e_EFFSF_DOWN = ElectronSF(EGEffCorrections,
                                                                                        event,
                                                                                        self.metadata.datatype,
                                                                                        pileup_tool,
                                                                                        YEAR, event.RunNumber, event.isLTT)

                #Lepton Trigger scale factors
                if not event.isLTT:
                    if event.leptonType == 'mu':
                        leptontrigsf_w, tree.sys_mu_TRIGSF_UP, tree.sys_mu_TRIGSF_DOWN = LeptonSLTSF(SLTCorrections,
                                                                                                     event,
                                                                                                     self.metadata.datatype,
                                                                                                     pileup_tool,
                                                                                                     YEAR)
                else:
                    if event.leptonType == 'mu':
                        leptontrigsf_w, tree.sys_mu_TRIGSF_UP, tree.sys_mu_TRIGSF_DOWN = MuonLTTSF(LTTCorrections,
                                                                                                   event,
                                                                                                   self.metadata.datatype,
                                                                                                   YEAR,
                                                                                                   pileup_tool)
                    if event.leptonType == 'e':
                        leptontrigsf_w, tree.sys_e_TRIGSF_UP, tree.sys_e_TRIGSF_DOWN = ElectronLTTSF(LTTCorrections,
                                                                                                     event,
                                                                                                     self.metadata.datatype,
                                                                                                     YEAR,
                                                                                                     pileup_tool,
                                                                                                     event.RunNumber)

                    #Tau trigger scale factors
                    tautriggersf_w, tree.sys_tau_TRIGSF_UP, tree.sys_tau_TRIGSF_DOWN = TauTriggerSF(LTTtauCorrections1,
                                                                                                    LTTtauCorrections2,
                                                                                                    LTTtauCorrections3,
                                                                                                    Tau,
                                                                                                    npileup_vtx,
                                                                                                    event.RunNumber,
                                                                                                    YEAR,
                                                                                                    event.leptonType,
                                                                                                    pileup_tool,
                                                                                                    self.metadata.datatype)
                if event.leptonType == 'mu':
                    muonisosf_w, tree.sys_mu_ISOSF_UP, tree.sys_mu_ISOSF_DOWN = MuonIsoEffCorrection(MIEffCorrections,
                                                                                                     event,
                                                                                                     self.metadata.datatype,
                                                                                                     YEAR)

                event_weight = mc_w * tautriggersf_w * tauidsf_w * leptonsf_w * leptontrigsf_w * muonisosf_w

                tree.w_mc           = mc_w
                tree.w_tauidsf      = tauidsf_w
                tree.w_leptonsf     = leptonsf_w
                tree.w_leptontrigsf = leptontrigsf_w
                tree.w_tautriggersf = tautriggersf_w
                tree.w_muonisosf    = muonisosf_w

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

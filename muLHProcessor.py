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
from higgstautau.filters import PriVertex, JetCleaning, LArError, LArHole
from higgstautau import mass
from higgstautau.trigger import update_trigger_config, get_trigger_config
from higgstautau.trigger import utils as triggerutils
from higgstautau.pileup import PileupReweighting

from goodruns import GRL
import subprocess
    
class muLHProcessor(ATLASStudent):
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

        OutputModel = RecoTauMuBlock + EventVariables + RecoMET

        onfilechange = []
        if self.metadata.datatype == datasets.DATA:
            merged_grl = GRL()

            def update_grl(student, grl, name, file):
                grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())

            onfilechange.append((update_grl, (self, merged_grl,)))

        # update the trigger config maps on every file change
        onfilechange.append((update_trigger_config, (trigger_config,)))

        merged_cutflow = Hist(6, 0, 6, name='cutflow', type='D')
        
        def update_cutflow(student, cutflow, name, file):
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
        tree = Tree(name=self.metadata.name + '_mulh', model=OutputModel)


        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']
        
        if self.metadata.datatype == datasets.MC:
            copied_variables += mc_triggers

        tree.set_buffer(chain.buffer, variables=copied_variables, create_branches=True, visible=False)

        chain.always_read(copied_variables)

        # set the event filters
        # passthrough for MC for trigger acceptance studies
        event_filters = EventFilterList([
            GRLFilter(self.grl, passthrough=self.metadata.datatype != datasets.DATA),
            PriVertex(),
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
        chain.define_collection(name="electrons", prefix="el_", size="el_n", mix=FourMomentum)
        chain.define_collection(name="taus", prefix="tau_", size="tau_n", mix=FourMomentum)
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        chain.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        # define tree objects
        tree.define_object(name='tau', prefix='tau_')
        tree.define_object(name='muon', prefix='muon_')

        if self.metadata.datatype == datasets.MC:
            # Initialize the pileup reweighting tool
            pileup_tool = PileupReweighting()
            pileup_tool.AddConfigFile('higgstautau/pileup/%s_defaults.prw.root' % self.metadata.category)
            pileup_tool.AddLumiCalcFile('grl/lumicalc/lephad/ilumicalc_histograms_None_178044-191933.root')
            # discard unrepresented data (with mu not simulated in MC)
            pileup_tool.SetUnrepresentedDataAction(1)
            pileup_tool.Initialize()


        # entering the main event loop...
        for event in chain:

            tree.reset()
            cutflow.reset()
            
            # Select tau with highest BDT score and surviving muon
            Tau = event.taus[0]
            Muon = event.muons[0]

            """
            RecoTauMuBlock filling
            """
            RecoTauMuBlock.set(event, tree, Tau, Muon)


            """
            Jets
            """
            numJets = len(event.jets)
            tree.numJets = numJets
            
            for jet in event.jets:
                tree.jet_fourvect.push_back(jet.fourvect)
                tree.jet_jvtxf.push_back(jet.jvtxf)


            """
            Miscellaneous
            """
            tree.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or
                                         (vtx.type == 3 and vtx.nTracks >= 2)])

            HT = 0

            HT += Tau.fourvect.Pt()
            HT += Muon.fourvect.Pt()
            for jet in event.jets:
                HT += jet.fourvect.Pt()

            tree.HT = HT

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
            muPhiVector = Vector2(Muon.fourvect.Px(), Muon.fourvect.Py())
            mT = sqrt(2*MET*muET*(1 - cos(muPhiVector.DeltaPhi(MET_vect))))
            tree.mass_transverse_met_muon = mT

            """
            Higgs fancier mass calculation
            """
            collin_mass, tau_x, muon_x = mass.collinearmass(Tau, Muon, METx, METy)
            tree.mass_collinear_tau_muon = collin_mass
            tree.tau_x  = tau_x
            tree.muon_x = muon_x
            tree.mass_mmc_tau_muon = 91#mass.missingmass(Tau, Muon, event.jets, METx, METy, sumET, self.metadata.datatype, 1)



            """
            Calculate fancier quantities for BDT input
            """

            mass_j1_j2 = -1111
            eta_product_j1_j2 = -1111
            eta_delta_j1_j2 = -1111
            tau_centrality_j1_j2 = -1111
            muon_centrality_j1_j2 = -1111
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

            PtSum  += (Tau.fourvect.Pt() + Muon.fourvect.Pt())
            PtSum2 += (Tau.fourvect.Pt()**2 + Muon.fourvect.Pt()**2)

            tree.neff_pt = PtSum2/(PtSum**2)
            tree.mass_all_jets = allJets.M()
            
            if len(event.jets) >= 2:
                event.jets.sort(key=lambda jet: jet.pt, reverse=True)
                
                jet1 = event.jets[0].fourvect
                jet2 = event.jets[1].fourvect

                mass_j1_j2 = (jet1 + jet2).M()
                eta_product_j1_j2 = jet1.Eta() * jet2.Eta()
                eta_delta_j1_j2 = abs(jet1.Eta() - jet2.Eta())
                tau_centrality_j1_j2 = eventshapes.eta_centrality(Tau.fourvect.Eta(), jet1.Eta(), jet2.Eta())
                muon_centrality_j1_j2 = eventshapes.eta_centrality(Muon.fourvect.Eta(), jet1.Eta(), jet2.Eta())

                sphericity, aplanarity = eventshapes.sphericity_aplanarity([Tau.fourvect,
                                                                            Muon.fourvect,
                                                                            jet1,
                                                                            jet2])


            tree.mass_j1_j2 = mass_j1_j2
            tree.eta_product_j1_j2 = eta_product_j1_j2
            tree.eta_delta_j1_j2 = eta_delta_j1_j2
            tree.tau_centrality_j1_j2 = tau_centrality_j1_j2
            tree.muon_centrality_j1_j2 = muon_centrality_j1_j2

            tau2Vector = Vector2(Tau.fourvect.Px(), Tau.fourvect.Py())
            muon2Vector = Vector2(Muon.fourvect.Px(), Muon.fourvect.Py())
            
            tree.met_phi_centrality = eventshapes.phi_centrality(tau2Vector, muon2Vector, MET_vect)

            tree.sphericity = sphericity
            tree.aplanarity = aplanarity

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
                    
            tree.weight = event_weight

            # fill output ntuple
            tree.cutflow = cutflow.int()
            tree.Fill()

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

        if self.metadata.datatype == datasets.DATA:
            xml_string = ROOT.TObjString(merged_grl.str())
            xml_string.Write('lumi')
        merged_cutflow.Write()

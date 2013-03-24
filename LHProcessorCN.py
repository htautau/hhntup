import ROOT, sys
import math, random

from argparse import ArgumentParser

from rootpy.tree.filtering import *
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.tree.cutflow import Cutflow
from rootpy.math.physics.vector import Vector2, LorentzVector
from rootpy.plotting import Hist
from rootpy.io import open as ropen

from higgstautau.mixins import *
from higgstautau.lephad.models import *
from higgstautau.filters import *
from higgstautau.lephad.filters import *
from higgstautau import eventshapes
from higgstautau.systematics import Systematics

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.batch import ATLASStudent

from goodruns import GRL
import subprocess

YEAR = 2011
VERBOSE = False

## Shape systematics dictionary to access systematics trees
SYSTEMATICS = {
    ## JES systematics
    'ATLAS_JES_BASE_DOWN' : 'SystematicsDOWN/JES_BASE',
    'ATLAS_JES_BASE_UP'   : 'SystematicsUP/JES_BASE',
    'ATLAS_JES_FLAV_DOWN' : 'SystematicsDOWN/JES_FLAV',
    'ATLAS_JES_FLAV_UP'   : 'SystematicsUP/JES_FLAV',
    'ATLAS_JES_FWD_DOWN'  : 'SystematicsDOWN/JES_FWD',
    'ATLAS_JES_FWD_UP'    : 'SystematicsUP/JES_FWD',

    ## TES systematics
    'ATLAS_TAU_ES_DOWN' : 'SystematicsDOWN/TES',
    'ATLAS_TAU_ES_UP'   : 'SystematicsUP/TES',

    ## MET systematics
    'ATLAS_MET_RESOSOFT_DOWN'  : 'SystematicsDOWN/METResSys',
    'ATLAS_MET_RESOSOFT_UP'    : 'SystematicsUP/METResSys',
    'ATLAS_MET_SCALESOFT_DOWN' : 'SystematicsDOWN/METScaleSys',
    'ATLAS_MET_SCALESOFT_UP'   : 'SystematicsUP/METScaleSys',

    ## Electron systematics
    'ATLAS_EL_ES_DOWN'  : 'SystematicsDOWN/EESSys',
    'ATLAS_EL_ES_UP'    : 'SystematicsUP/EESSys',
    'ATLAS_EL_RES_DOWN' : 'SystematicsDOWN/ElEnResSys',
    'ATLAS_EL_RES_UP'   : 'SystematicsUP/ElEnResSys',

    ## Muon systematics
    'ATLAS_MU_ES_DOWN' : 'SystematicsDOWN/MuSys',
    'ATLAS_MU_ES_UP'   : 'SystematicsUP/MuSys'
}

## Load overlap events
f = open('overlap.txt')
lines = f.readlines()

overlap = []

for line in lines:
    elements = line.split()
    if elements[0] == '<':
        overlap.append('.'.join(elements[1:]))

class LHProcessorCN(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def __init__(self, options, **kwargs):

        super(LHProcessorCN, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--syst-terms', default=None)
        self.args = parser.parse_args(options)
        if self.args.syst_terms is not None:
            self.args.syst_terms = self.args.syst_terms.split(',')

    @staticmethod
    def merge(inputs, output, metadata):

        # merge output trees
        root_output = output + '.root'
        subprocess.call(['hadd', root_output] + inputs)

        # if metadata.datatype == datasets.DATA:
        #     # merge GRLs
        #     grl = GRL()
        #     for input in inputs:
        #         grl |= GRL('%s:/lumi' % input)
        #     grl.save('%s:/lumi' % root_output)

            
    def work(self):
        """
        This is the one function that all "ATLASStudent"s must implement.
        """

        ## Use the metadata to set the year
        YEAR = self.metadata.year

        ## Set the output ntuple model
        OutputModel = RecoTauLepBlock + EventVariables + RecoMET# + SysWeights

        onfilechange = []
        #if self.metadata.datatype == datasets.DATA:
            #merged_grl = GRL()

            #def update_grl(student, grl, name, file, tree):
            #    grl |= str(file.Get('Lumi/%s' % student.metadata.treename).GetString())
            #    
            #onfilechange.append((update_grl, (self, merged_grl,)))

        if self.metadata.datatype != datasets.EMBED:
            merged_cutflow = Hist(1, 0, 1, name='TotalEvents', type='D')

            def update_cutflow(student, cutflow, name, file, tree):
                cutflow[0] += file.TotalEvents[0]
                
            onfilechange.append((update_cutflow, (self, merged_cutflow,)))

        ## initialize the TreeChain of all input files

        ## Nominal
        if self.args.syst_terms is None:
            chain = TreeChain(self.metadata.treename,
                              files=self.files,
                              events=self.events,
                              cache=True,
                              cache_size=10000000,
                              learn_entries=30,
                              onfilechange=onfilechange,
                              verbose=True)

        ## Shape systematic variation
        elif len(self.args.syst_terms) == 1:
            chain = TreeChain(SYSTEMATICS[self.args.syst_terms[0]],
                              files=self.files,
                              events=self.events,
                              cache=True,
                              cache_size=10000000,
                              learn_entries=30,
                              onfilechange=onfilechange,
                              verbose=True)
                
        else:
            print "ERROR: Too many systematics terms, LHProcessorCN can only handle one at a time."
            sys.exit()

        # create output tree
        self.output.cd()
        tree = Tree(name='lh', model=OutputModel)


        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing',
                            'RunNumber',
                            'EventNumber',
                            'lbn']

        tree.set_buffer(chain.buffer, branches=copied_variables, create_branches=True, visible=False)

        chain.always_read(copied_variables)

        ## Setting event filters
        event_filters = EventFilterList([
            JetSelection(YEAR, False),
            #VBFFilter(2.1, 210),
            #AntiVBFFilter(4.1, 410),
        ])

        self.filters['event'] = event_filters

        chain.filters += event_filters

        if self.metadata.datatype != datasets.EMBED:
            cutflow = Cutflow()

        ## Define handle to jet collection
        chain.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentumMeV)
        chain.define_collection(name="truthjets", prefix="jet_antikt4truth_", size="jet_antikt4truth_n", mix=FourMomentum)
            
        # define tree objects
        tree.define_object(name='tau', prefix='tau_')
        tree.define_object(name='lep', prefix='lep_')

        # entering the main event loop...
        for event in chain:

            if self.metadata.datatype != datasets.EMBED:
                cutflow.reset()

            ## Take care of overlap events for data
            if self.metadata.datatype == datasets.DATA:
                current_event = '%d.%d' % (event.RunNumber, event.EventNumber)
                if current_event in overlap:
                    overlap[overlap.index(current_event)] += '.done'
                    continue
                
            ## -- Subsample -- ##
                
            ## Initialize the subsample flags
            tree.subsample1 = False
            tree.subsample2 = False
            tree.subsample2_1 = False
            tree.subsample2_2 = False

            ## Determine the subsamples
            subsample1_1 = (event.EventNumber%4 == 0)
            subsample1_2 = (event.EventNumber%4 == 2)
            subsample2_1 = (event.EventNumber%4 == 1)
            subsample2_2 = (event.EventNumber%4 == 3)

            ## Assign the subsamples
            if subsample1_1 or subsample1_2: tree.subsample1 = True
            if subsample2_1 or subsample2_2: tree.subsample2 = True
            tree.subsample2_1 = subsample2_1
            tree.subsample2_2 = subsample2_2

            ## -- Event selection -- ##
            leptype = None
            if event.evtsel_is_mutau or event.evtsel_is_mu:
                leptype = 0
            elif event.evtsel_is_eltau or event.evtsel_is_el:
                leptype = 1
            else:
                continue

            #if not event.evtsel_is_tau: continue
            tree.is_tau = evtsel_is_tau
            if self.metadata.datatype != datasets.DATA:
                if not event.evtsel_is_tau: continue
            #if not event.evtsel_is_dilepVeto: continue

            ## Convert to private ntuple
            ##########################################################################

            ## Event type flags
            tree.dilep_veto = event.evtsel_is_dilepVeto
            tree.dilep_control = event.evtsel_is_dilepCR

            LTT = (event.evtsel_is_mutau or event.evtsel_is_eltau)
            SLT = (event.evtsel_is_mu or event.evtsel_is_el)

            tree.LTT = (LTT and not SLT)
            tree.SLT = (SLT and not LTT)

            
            ## -- Tau Information -- ##
            ## Basic kinematics
            Tau = LorentzVector()
            TauPt  = event.evtsel_tau_et*GeV
            TauEta = event.evtsel_tau_eta
            TauPhi = event.evtsel_tau_phi
            Tau.SetPtEtaPhiM(TauPt, TauEta, TauPhi, 0)
            getattr(tree, 'tau_fourvect').set_from(Tau)

            TrueTau = LorentzVector()
            TrueTauPt  = event.evtsel_truethad_pt
            TrueTauEta = event.evtsel_truethad_eta
            TrueTauPhi = event.evtsel_truethad_phi
            TrueTau.SetPtEtaPhiM(TrueTauPt, TrueTauEta, TrueTauPhi, 0)
            

            ## Other decorations
            setattr(tree, 'tau_BDTJetScore', event.evtsel_tau_JetBDTScore)
            setattr(tree, 'tau_BDTEleScore', -1)
            setattr(tree, 'tau_JetBDTSigLoose', 0)
            setattr(tree, 'tau_JetBDTSigMedium', event.evtsel_tau_is_Medium)
            setattr(tree, 'tau_JetBDTSigTight', event.evtsel_tau_is_Tight)
            setattr(tree, 'tau_EleBDTtight', event.evtsel_tau_eVetoTight)
            setattr(tree, 'tau_numTrack', event.evtsel_tau_numTrack)
            setattr(tree, 'tau_charge', event.evtsel_tau_charge)
            setattr(tree, 'tau_isTrueLep', event.evtsel_tau_is_el or event.evtsel_tau_is_mu)

            

            ## -- Lepton Information -- ##

            ## Basic kinematics
            Lep = LorentzVector()
            LepPt  = event.evtsel_lep_pt*GeV
            LepEta = event.evtsel_lep_eta
            LepPhi = event.evtsel_lep_phi
            Lep.SetPtEtaPhiM(LepPt, LepEta, LepPhi, 0)
            getattr(tree, 'lep_fourvect').set_from(Lep)

            TrueLep = LorentzVector()
            TrueLepPt  = event.evtsel_truetlep_pt
            TrueLepEta = event.evtsel_truetlep_eta
            TrueLepPhi = event.evtsel_truetlep_phi
            TrueLep.SetPtEtaPhiM(TrueLepPt, TrueLepEta, TrueLepPhi, 0)

            ## Other information
            setattr(tree, 'lep_charge', event.evtsel_lep_charge)
            setattr(tree, 'lep_isolated', event.evtsel_is_isoLep)
            setattr(tree, 'lep_leptype', leptype)


            
            ## -- Jet Information -- ##
            
            tree.numJets = len(event.jets)

            numJets50 = 0
            numJets30 = 0

            tree.btag = False

            vector_all = Tau + Lep
            
            for jet in event.jets:
                vector_all += jet.fourvect
                #tree.jet_fourvect.push_back(jet.fourvect)
                #tree.jet_jvtxf.push_back(jet.jvtxf)
                #tree.jet_btag.push_back(jet.flavor_weight_MV1)
                if jet.flavor_weight_MV1 > 0.7892:
                    tree.btag = True
                if jet.fourvect.Pt() > 30*GeV:
                    numJets30 += 1
                    if jet.fourvect.Pt() > 50*GeV:
                        numJets50 += 1

            tree.leadjet_btag     = -1111
            tree.leadjet_jvtxf    = -1111

            tree.subleadjet_btag  = -1111
            tree.subleadjet_jvtxf = -1111
                        
            if tree.numJets > 0:
                tree.leadjet_fourvect = event.jets[0].fourvect
                tree.leadjet_btag     = event.jets[0].flavor_weight_MV1
                tree.leadjet_jvtxf    = event.jets[0].jvtxf

                if tree.numJets > 1:
                    tree.subleadjet_fourvect = event.jets[1].fourvect
                    tree.subleadjet_btag     = event.jets[1].flavor_weight_MV1
                    tree.subleadjet_jvtxf    = event.jets[1].jvtxf

                    

            tree.numJets50 = numJets50
            tree.numJets30 = numJets30

            # ## true jets
            # if ( self.metadata.datatype == datasets.MC ):
            #     for truthjet in event.truthjets:
            #         tree.truthjet_fourvect.push_back(truthjet.fourvect)
                    

            
            ## -- MET Information -- ##

            METx = event.evtsel_met_etx*GeV
            METy = event.evtsel_met_ety*GeV
            MET_vect = Vector2(METx, METy)
            MET = MET_vect.Mod()
            tree.MET = MET
            getattr(tree, 'MET_vect').set_from(MET_vect)

            tree.MET_Reftau_pt = event.evtsel_MET_RefTau_sumet*GeV


            
            ## -- Variables -- ##

            #transverse mass
            lepPhiVector = Vector2(Lep.Px(), Lep.Py())
            tauPhiVector = Vector2(Tau.Px(), Tau.Py())
            dPhi_MET_lep = lepPhiVector.DeltaPhi(MET_vect)
            dPhi_MET_tau  = tauPhiVector.DeltaPhi(MET_vect)
            mT = sqrt(2*MET*LepPt*(1 - cos(dPhi_MET_lep)))
            mTtau = sqrt(2*MET*TauPt*(1 - cos(dPhi_MET_tau)))
            tree.mass_transverse_met_lep = mT
            tree.mass_transverse_met_tau = mTtau
            tree.dphi_met_lep = dPhi_MET_lep

            tree.pt_vector_sum_all = vector_all.Pt()

            #ddR
            tree.ddr_tau_lep, tree.dr_tau_lep, tree.resonance_pt_tau_lep = eventshapes.DeltaDeltaR(Tau, Lep, MET_vect)

            ## Weird
            tree.MET_over_resonance_pt = MET/(tree.resonance_pt_tau_lep*GeV)
            tree.sum_dphi_tau_lep = event.evtsel_dPhiSum
            tree.dpt_tau_lep = Lep.Pt() - Tau.Pt()
            
            #number of vertices
            tree.nvtx = event.evtsel_vertices

            #sumPt
            sumPt = 0
            sumJetPt = 0
            sumPt += TauPt
            sumPt += LepPt
            for jet in event.jets:
                sumPt += jet.fourvect.Pt()
                sumJetPt += jet.fourvect.Pt()

            tree.sumPt = sumPt
            tree.HT_jets = sumJetPt

            #Ditau quantities
            tree.tau_x  = event.evtsel_tau_x1
            tree.lep_x = event.evtsel_lep_x2
            tau_x_lep_x = event.evtsel_tau_x1*event.evtsel_lep_x2
            tree.mass_mmc_tau_lep = event.evtsel_ditau_MMC
            tree.pt_mmc_tau_lep = -1
            tree.met_mmc_tau_lep = -1
            tree.mass_vis_tau_lep = event.evtsel_ditau_visibleMass*GeV
            tree.charge_product_tau_lep = event.evtsel_tau_charge * event.evtsel_lep_charge
            tree.pt_ratio_tau_lep = LepPt/TauPt
            tree.met_phi_centrality = eventshapes.phi_centrality(tauPhiVector, lepPhiVector, MET_vect)
            tree.pt_balance_tau_lep = (LepPt-TauPt)/(LepPt+TauPt)

            #VBF topology
            mass_j1_j2 = -1111
            eta_product_j1_j2 = -1111
            eta_delta_j1_j2 = -1111
            eta_balance_j1_j2 = -1111
            tau_centrality_j1_j2 = -1111
            lep_centrality_j1_j2 = -1111
            tau_lep_centrality_j1_j2 = -1111
            sphericity = -1111
            min_deta_tau_lep_j1_j2 = -1111
            true_dphi_resonance_dijet = -1111

            event.jets.sort(key=lambda jet: jet.fourvect.Pt(), reverse=True)
                
            leadJetPt = 0.0
            if len(event.jets) > 0:
                leadJetPt = event.jets[0].fourvect.Pt()
            tree.leadJetPt = leadJetPt

            if len(event.jets) >= 2:

                jet1Pt  = event.jets[0].fourvect.Pt()
                jet1Eta = event.jets[0].fourvect.Eta()
                jet1Phi = event.jets[0].fourvect.Phi()
                jet1 = LorentzVector()
                jet1.SetPtEtaPhiM(jet1Pt, jet1Eta, jet1Phi, 0.0)

                ## Truth match
                trueJet1 = LorentzVector()
                for truthjet in event.truthjets:
                    if jet1.DeltaR(truthjet.fourvect) < 0.2:
                        trueJet1 = truthjet.fourvect
                        

                jet2Pt  = event.jets[1].fourvect.Pt()
                jet2Eta = event.jets[1].fourvect.Eta()
                jet2Phi = event.jets[1].fourvect.Phi()
                jet2 = LorentzVector()
                jet2.SetPtEtaPhiM(jet2Pt, jet2Eta, jet2Phi, 0.0)

                ## Truth match
                trueJet2 = LorentzVector()
                for truthjet in event.truthjets:
                    if jet2.DeltaR(truthjet.fourvect) < 0.2:
                        trueJet2 = truthjet.fourvect

                jet1_2Vector = Vector2(jet1.Px(), jet1.Py())
                jet2_2Vector = Vector2(jet2.Px(), jet2.Py())

                mass_j1_j2 = (jet1 + jet2).M()
                eta_product_j1_j2 = jet1.Eta() * jet2.Eta()
                eta_delta_j1_j2 = abs(jet1.Eta() - jet2.Eta())
                eta_balance_j1_j2 = (jet1.Eta()**2 + jet2.Eta()**2)*abs(jet2.Eta() - jet1.Eta())
                
                tau_centrality_j1_j2 = eventshapes.eta_centrality(Tau.Eta(), jet1.Eta(), jet2.Eta())
                lep_centrality_j1_j2 = eventshapes.eta_centrality(Lep.Eta(), jet1.Eta(), jet2.Eta())
                tau_lep_centrality_j1_j2 = tau_centrality_j1_j2 * lep_centrality_j1_j2

                deta_tau_lep_j1 = abs( (Tau + Lep).Eta() - jet1.Eta() )
                deta_tau_lep_j2 = abs( (Tau + Lep).Eta() - jet2.Eta() )
                
                min_deta_tau_lep_j1_j2 = min(deta_tau_lep_j1, deta_tau_lep_j2)

                ## True information
                higgs = (TrueTau + TrueLep)
                higgs_2Vector = Vector2(higgs.Px(), higgs.Py())
                dijet = (trueJet1 + trueJet2)
                dijet_2Vector = Vector2(dijet.Px(), dijet.Py())

                true_dphi_resonance_dijet = higgs_2Vector.DeltaPhi(dijet_2Vector)

            jets_vector_sum = LorentzVector()

            allJetList = []
            for jet in event.jets:
                JetPt  = jet.fourvect.Pt()
                JetEta = jet.fourvect.Eta()
                JetPhi = jet.fourvect.Phi()
                Jet = LorentzVector()
                Jet.SetPtEtaPhiM(JetPt, JetEta, JetPhi, 0.0)
                jets_vector_sum += Jet
                allJetList.append(Jet)
            
            sphericity, aplanarity = eventshapes.sphericity_aplanarity([Tau, Lep] + allJetList)

            tree.mass_j1_j2 = mass_j1_j2
            
            tree.eta_product_j1_j2 = eta_product_j1_j2
            tree.eta_delta_j1_j2 = eta_delta_j1_j2
            tree.eta_balance_j1_j2 = eta_balance_j1_j2
            
            tree.tau_centrality_j1_j2 = tau_centrality_j1_j2
            tree.lep_centrality_j1_j2 = lep_centrality_j1_j2
            tree.tau_lep_centrality_j1_j2 = tau_lep_centrality_j1_j2

            tree.sphericity = sphericity

            tree.min_deta_tau_lep_j1_j2 = min_deta_tau_lep_j1_j2
            tree.true_dphi_resonance_dijet = true_dphi_resonance_dijet



            ## -- Categories -- ##
            tree.category_vbf_train = False
            tree.category_vbf_test = False
            tree.category_boosted = False
            tree.category_1j = False
            tree.category_0j = False

            if (numJets30 >= 2 and numJets50 >= 1):
                tree.category_vbf_train = True
            if (numJets30 >= 2 and numJets50 >= 1 and eta_delta_j1_j2 > 3.0):
                tree.category_vbf_test = True
            elif tree.resonance_pt_tau_lep > 100:
                tree.category_boosted = True
            elif tree.numJets >= 1:
                tree.category_1j = True
            elif tree.numJets == 0:
                tree.category_0j = True



            ## -- Event weights -- ##
            tree.weight = event.evtsel_weight


            ## -- True Information -- ##
            tree.true_higgs_mass = (TrueTau + TrueLep).M()
            
            

            # fill output ntuple
            if self.metadata.datatype != datasets.EMBED:
                tree.cutflow = cutflow.int()
            tree.Fill(reset=True)

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

        #if self.metadata.datatype == datasets.DATA:
        #    xml_string = ROOT.TObjString(merged_grl.str())
        #    xml_string.Write('lumi')
        #if self.metadata.datatype == datasets.MC:
        #    merged_cutflow.Write()

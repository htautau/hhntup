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
DO_VBF_OLR = False
VBF_OLR_CONFIG = {
    'ZmumuNp2' : (False, -1, -1, True, 2.1, 210),
    'ZmumuNp3' : (False, -1, -1, True, 2.1, 210),
    'ZmumuNp4' : (False, -1, -1, True, 2.1, 210),
    'ZmumuNp5' : (False, -1, -1, True, 2.1, 210),
    'ZeeNp2' : (False, -1, -1, True, 2.1, 210),
    'ZeeNp3' : (False, -1, -1, True, 2.1, 210),
    'ZeeNp4' : (False, -1, -1, True, 2.1, 210),
    'ZeeNp5' : (False, -1, -1, True, 2.1, 210),
    'ZtautauNp2' : (False, -1, -1, True, 2.1, 210),
    'ZtautauNp3' : (False, -1, -1, True, 2.1, 210),
    'ZtautauNp4' : (False, -1, -1, True, 2.1, 210),
    'ZtautauNp5' : (False, -1, -1, True, 2.1, 210),
    
    'WtaunuNp3' : (False, -1, -1, True, 3.1, 410),
    'WtaunuNp4' : (False, -1, -1, True, 3.1, 410),
    'WtaunuNp5' : (False, -1, -1, True, 3.1, 410),
    'WenuNp3' : (False, -1, -1, True, 3.1, 410),
    'WenuNp4' : (False, -1, -1, True, 3.1, 410),
    'WenuNp5' : (False, -1, -1, True, 3.1, 410),
    'WmunuNp3' : (False, -1, -1, True, 3.1, 410),
    'WmunuNp4' : (False, -1, -1, True, 3.1, 410),
    'WmunuNp5' : (False, -1, -1, True, 3.1, 410),

    
    'VBF_ZmumuNp2' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZmumuNp3' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZmumuNp4' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZmumuNp5' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZeeNp2' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZeeNp3' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZeeNp4' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ZeeNp5' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ATau_ZtautauNp2' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ATau_ZtautauNp3' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ATau_ZtautauNp4' : (True, 2.1, 210, True, 4.1, 410),
    'VBF_ATau_ZtautauNp5' : (True, 2.1, 210, True, 4.1, 410),

    'VBF_WtaunuNp3' : (True, 3.1, 410, False, -1, -1),
    'VBF_WtaunuNp4' : (True, 3.1, 410, False, -1, -1),
    'VBF_WtaunuNp5' : (True, 3.1, 410, False, -1, -1),
    'VBF_WenuNp3' : (True, 3.1, 410, False, -1, -1),
    'VBF_WenuNp4' : (True, 3.1, 410, False, -1, -1),
    'VBF_WenuNp5' : (True, 3.1, 410, False, -1, -1),
    'VBF_WmunuNp3' : (True, 3.1, 410, False, -1, -1),
    'VBF_WmunuNp4' : (True, 3.1, 410, False, -1, -1),
    'VBF_WmunuNp5' : (True, 3.1, 410, False, -1, -1),

    'TVBF_ZmumuNp2' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZmumuNp3' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZmumuNp4' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZmumuNp5' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZeeNp2' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZeeNp3' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZeeNp4' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ZeeNp5' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ATau_ZtautauNp2' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ATau_ZtautauNp3' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ATau_ZtautauNp4' : (True, 4.1, 410, False, -1, -1),
    'TVBF_ATau_ZtautauNp5' : (True, 4.1, 410, False, -1, -1),
}

## Shape systematics dictionary to access systematics trees
SYSTEMATICS = {
    ## JES systematics
    'ATLAS_JES_BASE_DOWN' : 'SystematicsDOWN/JES_BASE',
    'ATLAS_JES_BASE_UP'   : 'SystematicsUP/JES_BASE',
    'ATLAS_JES_FLAV_DOWN' : 'SystematicsDOWN/JES_FLAV',
    'ATLAS_JES_FLAV_UP'   : 'SystematicsUP/JES_FLAV',
    'ATLAS_JES_FWD_DOWN'  : 'SystematicsDOWN/JES_FWD',
    'ATLAS_JES_FWD_UP'    : 'SystematicsUP/JES_FWD',
    'ATLAS_JER_DOWN'      : 'SystematicsDOWN/JER',
    'ATLAS_JER_UP'        : 'SystematicsUP/JER',

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

        filter_list = [JetSelection(YEAR)]


        if DO_VBF_OLR:
            vbf_filter_config = None

            try:
                vbf_filter_config = VBF_OLR_CONFIG[self.metadata.name]
            except KeyError:
                pass

            if vbf_filter_config is not None:
                if vbf_filter_config[0]:
                    filter_list.append(VBFFilter(vbf_filter_config[1], vbf_filter_config[2]))
                if vbf_filter_config[3]:
                    filter_list.append(VBFFilter(vbf_filter_config[4], vbf_filter_config[5]))
            

        ## Setting event filters
        event_filters = EventFilterList(filter_list)

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

            tree.is_tau = event.evtsel_is_tau
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
            
            tree.numJets = event.evtsel_jets_num

            numJets50 = 0
            numJets30 = 0

            tree.btag = event.evtsel_MV1

            vector_all = LorentzVector()

            sorted_jets = []
            
            for jet in event.jets:
                sorted_jets.append(jet.fourvect)
                if jet.fourvect.Pt() > 30*GeV:
                    numJets30 += 1
                    if jet.fourvect.Pt() > 50*GeV:
                        numJets50 += 1

            sorted_jets.sort(key=lambda jet : jet.Pt(), reverse=True)
            lead    = sorted_jets[0]
            sublead = sorted_jets[1]

            # if tree.numJets > 0:
            #     lead = LorentzVector()
            #     lead.SetPtEtaPhiM(event.evtsel_jet_leading_pt*GeV,
            #                       event.evtsel_jet_leading_eta,
            #                       event.evtsel_jet_leading_phi,
            #                       event.evtsel_jet_leading_m*GeV)
            #     tree.leadjet_fourvect = lead

            #     if tree.numJets > 1:
            #         sublead = LorentzVector()
            #         sublead.SetPtEtaPhiM(event.evtsel_jet_subleading_pt*GeV,
            #                              event.evtsel_jet_subleading_eta,
            #                              event.evtsel_jet_subleading_phi,
            #                              event.evtsel_jet_subleading_m*GeV)
            #         tree.subleadjet_fourvect = sublead
            #         vector_all += Tau + Lep + lead + sublead

            

            tree.numJets50 = numJets50
            tree.numJets30 = numJets30

            ## -- MET Information -- ##

            METx = event.evtsel_met_etx*GeV
            METy = event.evtsel_met_ety*GeV
            MET_vect = Vector2(METx, METy)
            MET = MET_vect.Mod()
            tree.MET = MET
            getattr(tree, 'MET_vect').set_from(MET_vect)
            
            if tree.numJets > 1:
                MET_4vect = LorentzVector()
                MET_4vect.SetPtEtaPhiM(MET_vect.Mod(), 0, MET_vect.Phi(), 0)
                vector_all += MET_4vect

            
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

            #Higgs Pt
            tree.ddr_tau_lep, tree.dr_tau_lep, tree.resonance_pt_tau_lep = eventshapes.DeltaDeltaR(Tau, Lep, MET_vect)
            

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
                
            if len(event.jets) > 0:
                tree.leadJetPt = lead.Pt()

            if len(event.jets) >= 2:

                ## Truth match leading jet
                trueJet1 = LorentzVector()
                for truthjet in event.truthjets:
                    if lead.DeltaR(truthjet.fourvect) < 0.2:
                        trueJet1 = truthjet.fourvect
                        

                ## Truth match subleading jet
                trueJet2 = LorentzVector()
                for truthjet in event.truthjets:
                    if sublead.DeltaR(truthjet.fourvect) < 0.2:
                        trueJet2 = truthjet.fourvect

                jet1_2Vector = Vector2(lead.Px(), lead.Py())
                jet2_2Vector = Vector2(sublead.Px(), sublead.Py())

                mass_j1_j2 = (lead + sublead).M()
                eta_product_j1_j2 = lead.Eta() * sublead.Eta()
                eta_delta_j1_j2 = abs(lead.Eta() - sublead.Eta())
                
                tau_centrality_j1_j2 = eventshapes.eta_centrality(Tau.Eta(), lead.Eta(), sublead.Eta())
                lep_centrality_j1_j2 = eventshapes.eta_centrality(Lep.Eta(), lead.Eta(), sublead.Eta())
                tau_lep_centrality_j1_j2 = tau_centrality_j1_j2 * lep_centrality_j1_j2


                ## True information
                higgs = (TrueTau + TrueLep)
                higgs_2Vector = Vector2(higgs.Px(), higgs.Py())
                dijet = (trueJet1 + trueJet2)
                dijet_2Vector = Vector2(dijet.Px(), dijet.Py())

                true_dphi_resonance_dijet = higgs_2Vector.DeltaPhi(dijet_2Vector)


            ## Push dijet variables into the tree (default values if previous if statement was skipped)
            tree.mass_j1_j2                = mass_j1_j2
            tree.eta_product_j1_j2         = eta_product_j1_j2
            tree.eta_delta_j1_j2           = eta_delta_j1_j2
            tree.tau_centrality_j1_j2      = tau_centrality_j1_j2
            tree.lep_centrality_j1_j2      = lep_centrality_j1_j2
            tree.tau_lep_centrality_j1_j2  = tau_lep_centrality_j1_j2
            tree.true_dphi_resonance_dijet = true_dphi_resonance_dijet


            if abs(eta_delta_j1_j2 - event.evtsel_jets_deltaEta) > 0.01 and eta_delta_j1_j2 > 0:
                print '-------------------------------------------------'
                print 'deltaEtajj mine = %.3f, common %.3f' % (eta_delta_j1_j2, event.evtsel_jets_deltaEta)

            if abs(eta_product_j1_j2 - event.evtsel_jetEta1_jetEta2) > 0.01 and eta_product_j1_j2 > 0:
                print 'productEtajj mine = %.3f, common %.3f' % (eta_product_j1_j2, event.evtsel_jetEta1_jetEta2)

            if abs(mass_j1_j2/GeV - event.evtsel_dijet_mass) > 0.01 and mass_j1_j2 > 0:
                print 'Mjj mine = %.3f, common %.3f' % (mass_j1_j2/GeV, event.evtsel_dijet_mass)
            

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

            if tree.category_vbf_test or tree.category_boosted:
                tree.weight *= event.evtsel_bjet_weight


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
        if self.metadata.datatype == datasets.MC:
            merged_cutflow.Write()

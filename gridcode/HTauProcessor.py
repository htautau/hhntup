import ROOT
import math
from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from rootpy.tree.filtering import *
from atlastools.filtering import GRLFilter
from filters import *
from atlastools.batch import ATLASStudent
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.tree.cutflow import Cutflow
from mixins import MCParticle, FourMomentum, TauFourMomentum
import hepmc
import tautools
from models import *
import missingmass

ROOT.gErrorIgnoreLevel = ROOT.kFatal


class HTauProcessor(ATLASStudent):
    """
    ATLASStudent inherits from rootpy.batch.Student.
    """

    def work(self):
        """
        This is the one function that all "ATLASStudent"s must implement.
        """
         
        D4PD_model = RecoTauBlock + RecoJetBlock + EventVariables
        
        if self.fileset.datatype == datasets.MC:
            # this tree will contain info pertaining to true tau decays
            # for possible use in the optimization of a missing mass calculator
            mc_tree = Tree(name = "tau_mc", model=TrueTau_MCBlock)
            
            # only create truth branches for MC
            D4PD_model += TrueTauBlock
            
            # add branches for VBF Higgs associated partons
            if self.fileset.name.startswith("VBFH"):
                D4PD_model += PartonBlock
        
        
        # initialize the TreeChain of all input files (each containing one tree named self.fileset.treename)
        tree = TreeChain(self.fileset.treename,
                         files=self.fileset.files,
                         events=self.events,
                         usecache=True,
                         cache_size=10000000,
                         learn_entries=30)
        
        # create output tree
        self.output.cd()
        D4PD = Tree(name=self.fileset.name, model=D4PD_model)
        
        if self.fileset.datatype == datasets.MC:
            
            # store triggers for MC
            copied_variables = tree.glob("EF_*") + \
                               tree.glob("L1_*") + \
                               tree.glob("L2_*")
            # do a verbatim copy of these branches from the input tree into the output tree
            copied_variables += tree.glob("jet_AntiKt4TopoEM_*")
            D4PD.set_buffer(tree.buffer, variables=copied_variables, create_branches=True, visible=False)
            tree.always_read(copied_variables)
        
        # set the event filters
        # passthrough for MC for trigger acceptance studies
        self.event_filters = EventFilterList([
            GRLFilter(self.grl, passthrough = self.fileset.datatype != datasets.DATA),
            Trigger(passthrough = self.fileset.datatype == datasets.MC),
            PriVertex(),
            LArError(),
            JetCleaningLoose(passthrough = self.fileset.datatype != datasets.DATA),
            JetCleaningMedium(passthrough = self.fileset.datatype != datasets.DATA),
            LArHole(),
            JetCrackVeto(),
            ElectronVeto(),
            MuonVeto()
        ])
        tree.filters += self.event_filters

        twogoodtaus = EventFilter(name='TwoGoodTaus')
        twogoodjets = EventFilter(name='TwoGoodJets')

        self.event_filters += [
            twogoodtaus,
            twogoodjets
        ]

        cutflow = Cutflow()

        # define tree collections
        tree.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        # jet_eta etc is AntiKt4LCTopo in tau-perf D3PDs
        tree.define_collection(name="jets", prefix="jet_AntiKt4TopoEM_", size="jet_AntiKt4TopoEM_n", mix=FourMomentum)
        tree.define_collection(name="truetaus", prefix="trueTau_", size="trueTau_n")
        tree.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        tree.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n")
        tree.define_collection(name="electrons", prefix="el_", size="el_n")
        tree.define_collection(name="vertices", prefix="vxp_", size="vxp_n")
         
        # entering the main event loop...
        for event in tree:
             
            """
            Need to get all MC tau final states to build ntuple for missing mass calculator pdfs
            """ 
            if self.fileset.datatype == datasets.MC:
                tau_decays = tautools.get_tau_decays(event)
                for decay in tau_decays:
                    hadronic = decay.hadronic
                    if hadronic:
                        mc_tree.hadronic = True
                        mc_tree.nprong = decay.nprong
                        mc_tree.npi0 = decay.npi0
                        mc_tree.nneutrals = decay.nneutrals
                        mc_tree.fourvect.set_from(decay.fourvect)
                        mc_tree.fourvect_vis.set_from(decay.fourvect_visible)
                        mc_tree.fourvect_miss.set_from(decay.fourvect_missing)
                        mc_tree.dR_tau_nu = decay.dR_tau_nu
                        mc_tree.dTheta3d_tau_nu = decay.dTheta3d_tau_nu
                        mc_tree.Fill(reset=True)
                
            """
            Tau selection
            """
            # Only consider taus with at least a calo seed and which have at least one track
            taus = [tau for tau in event.taus if tau.author != 2 and tau.seedCalo_numTrack > 0]
            # kinematic region
            taus = [tau for tau in taus if tau.pt > 20*GeV]
            # muon veto
            taus = [tau for tau in taus if tau.muonVeto == 0]
            # charge requirement
            taus = [tau for tau in taus if abs(tau.charge) == 1]
            # Did not reconstruct two candidates so skip event
            if len(taus) < 2:
                twogoodtaus.failed()
                continue
            twogoodtaus.passed()

            # Sort the taus by BDT score
            taus = sorted(taus, key=lambda tau: tau.BDTJetScore, reverse=True)
            # Take the two taus with the highest BDT score
            taus = taus[:2]
            
            """
            Jet selection
            """
            # kinematic region
            jets = [jet for jet in event.jets if jet.pt > 20*GeV]
            jets = [jet for jet in event.jets if abs(jet.emscale_eta) < 4.5]
            
            """
            Overlap removal between taus and jets
            """
            otherjets = []
            for jet in jets:
                matched = False
                for tau in taus:
                    # tau.jet_emscale_eta/phi not in SMWZ D3PDs... using seedCalo_eta/phi instead
                    if utils.dR(jet.emscale_eta, jet.emscale_phi, tau.seedCalo_eta, tau.seedCalo_phi) < .2:
                        matched = True
                        break
                if not matched:
                    otherjets.append(jet)
            jets = otherjets
            
            """
            Will filter jets by jet vertex fraction (JVF) here...
            """

            """
            VBF jet selection
            two highest pT jets sorted by eta
            """
            best_jets = sorted(sorted(jets, key=lambda jet: jet.pt, reverse=True)[:2], key=lambda jet: jet.eta)
            if len(best_jets) < 2:
                # if there are fewer than 2 other jets then skip event
                twogoodjets.failed()
                continue
            twogoodjets.passed()

            """
            MET
            """
            MET_LocHadTopo_etx = event.MET_LocHadTopo_etx_CentralReg + event.MET_LocHadTopo_etx_EndcapRegion + event.MET_LocHadTopo_etx_ForwardReg
            MET_LocHadTopo_ety = event.MET_LocHadTopo_ety_CentralReg + event.MET_LocHadTopo_ety_EndcapRegion + event.MET_LocHadTopo_ety_ForwardReg
            
            try:
                MET_MuonBoy_etx = event.MET_MuonBoy_etx_CentralReg + event.MET_MuonBoy_etx_EndcapRegion + event.MET_MuonBoy_etx_ForwardReg
                MET_MuonBoy_ety = event.MET_MuonBoy_ety_CentralReg + event.MET_MuonBoy_ety_EndcapRegion + event.MET_MuonBoy_ety_ForwardReg
            except AttributeError: # Above are missing in tauPerf D3PDs
                MET_MuonBoy_etx = event.MET_MuonBoy_etx
                MET_MuonBoy_ety = event.MET_MuonBoy_ety

            try:
                MET_RefMuon_Track_etx = event.MET_RefMuon_Track_etx_CentralReg + event.MET_RefMuon_Track_etx_EndcapRegion + event.MET_RefMuon_Track_etx_ForwardReg
                MET_RefMuon_Track_ety = event.MET_RefMuon_Track_ety_CentralReg + event.MET_RefMuon_Track_ety_EndcapRegion + event.MET_RefMuon_Track_ety_ForwardReg
            except AttributeError: # Above are missing in tauPerf D3PDs
                MET_RefMuon_Track_etx = event.MET_RefMuon_Track_etx
                MET_RefMuon_Track_ety = event.MET_RefMuon_Track_ety
            
            METx = MET_LocHadTopo_etx + MET_MuonBoy_etx - MET_RefMuon_Track_etx
            METy = MET_LocHadTopo_ety + MET_MuonBoy_ety - MET_RefMuon_Track_ety

            MET = math.sqrt(METx**2 + METy**2)
            
            D4PD.MET = MET
            if MET > 0:
                phi = math.asin(METy / MET)
            else:
                phi = -1111.
            D4PD.MET_phi = phi
            
            # HT TODO: Fix
            sumET = event.MET_LocHadTopo_sumet + event.MET_MuonBoy_sumet - event.MET_RefMuon_Track_sumet
            D4PD.HT = sumET

            """
            MMC and misc variables
            """
            D4PD.MMC_mass = missingmass.mass(taus, jets, METx, METy, sumET, self.fileset.datatype)
            D4PD.Mvis_tau1_tau2 = utils.Mvis(taus[0].Et, taus[0].seedCalo_phi, taus[1].Et, taus[1].seedCalo_phi)
            D4PD.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or (vtx.type == 3 and vtx.nTracks >= 2)])
            D4PD.numJets = len(jets)
            if self.fileset.datatype == datasets.MC:
                D4PD.mu = event.lbn

            """
            Jet variables
            """
            RecoJetBlock.set(D4PD, best_jets[0], best_jets[1])

            """
            Experimenting here....
            Need to match jets to VBF jets
            """ 
            if self.fileset.datatype == datasets.MC:
                if self.fileset.name.startswith("VBFH"):
                    # get partons (already sorted by eta in hepmc)
                    parton1, parton2 = hepmc.get_VBF_partons(event)
                    PartonBlock.set(D4PD, parton1, parton2)
                    for jet in event.jets:
                        if jet in jets:
                            D4PD.jet_AntiKt4TopoEM_matched_dR.push_back(
                                min(
                                    utils.dR(jet.eta, jet.phi, parton1.eta, parton1.phi),
                                    utils.dR(jet.eta, jet.phi, parton2.eta, parton2.phi)
                                    )
                                )
                        D4PD.jet_AntiKt4TopoEM_matched_dR.push_back(1111)
            
            """
            Reco tau variables
            """
            RecoTauBlock.set(D4PD, taus[0], taus[1])
            
            """
            Truth-matching
            presently not possible in SMWZ D3PDs
            """
            if self.fileset.datatype == datasets.MC and hasattr(event, "trueTau_n"):
                unmatched_reco = [1, 2]
                unmatched_truth = range(1, min(len(event.truetaus)+1, 3))
                if len(event.truetaus) > 2:
                    print "ERROR: too many true taus: %i" % len(event.truetaus)
                    for tau in event.truetaus:
                        print "truth (pT: %.4f, eta: %.4f, phi: %.4f)" % (tau.pt, tau.eta, tau.phi),
                        if tau.tauAssoc_index >= 0:
                            matched_tau = event.taus[tau.tauAssoc_index]
                            print " ==> reco (pT: %.4f, eta: %.4f, phi: %.4f)" % (matched_tau.pt, matched_tau.seedCalo_eta, matched_tau.seedCalo_phi),
                            print "dR = %.4f" % tau.tauAssoc_dr
                        else:
                            print ""
                    D4PD.error = True
                matched_truth = []
                for i, tau in zip((1, 2), taus):
                    matching_truth_index = tau.trueTauAssoc_index
                    if matching_truth_index >= 0:
                        unmatched_reco.remove(i)
                        # check that this tau / true tau was not previously matched
                        if matching_truth_index+1 not in unmatched_truth or \
                           matching_truth_index in matched_truth:
                            print "ERROR: match collision!"
                            D4PD.tau1_matched_collision = True
                            D4PD.tau2_matched_collision = True
                            D4PD.trueTau1_matched_collision = True
                            D4PD.trueTau2_matched_collision = True
                            D4PD.error = True
                        else:
                            unmatched_truth.remove(matching_truth_index+1)
                            matched_truth.append(matching_truth_index)
                            setattr(D4PD, "tau%i_matched" % i, 1)
                            setattr(D4PD, "tau%i_matched_dR" % i, tau.trueTauAssoc_dr)
                            setattr(D4PD, "trueTau%i_matched" % i, 1)
                            setattr(D4PD, "trueTau%i_matched_dR" % i, event.truetaus[matching_truth_index].tauAssoc_dr)
                            TrueTauBlock.set(D4PD, i, event.truetaus[matching_truth_index])
                
                for i, j in zip(unmatched_reco, unmatched_truth):
                    TrueTauBlock.set(D4PD, i, event.truetaus[j-1])
             
            # fill output ntuple
            # use reset=True to reset all variables to their defaults after the fill
            # to avoid any values from this event carrying over into the next
            D4PD.cutflow = cutflow.int()
            D4PD.Fill(reset=True)
            cutflow.reset()

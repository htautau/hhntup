import ROOT
import math
from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from rootpy.tree.filtering import *
from atlastools.filtering import GRLFilter
from filters import *
from atlastools.batch import ATLASStudent
from rootpy.tree import Tree, TreeBuffer, TreeChain, TreeModel
from rootpy.types import *
from mixins import MCParticle
import hepmc
import tautools

ROOT.gSystem.Load("libMissingMassCalculator.so")
from rootpy.utils.classfactory import generate
generate("vector<TLorentzVector>", "<vector>;TLorentzVector.h")

ROOT.gErrorIgnoreLevel = ROOT.kFatal

class TruthTau(TreeModel):

    hadronic = Bool(default=False)
    nprong = Int(default=-1111)
    npi0 = Int(default=-1111)
    nneutrals = Int(default=-1111)
    pt = Float(default=0)
    eta = Float(default=-1111)
    phi = Float(default=-1111)
    pt_vis = Float(default=0)
    eta_vis = Float(default=-1111)
    phi_vis = Float(default=-1111)
    nu_pt = Float(default=0)
    nu_eta = Float(default=-1111)
    nu_phi = Float(default=-1111)
    dR_tau_nu = Float(default=-1111)
    dTheta3d_tau_nu = Float(default=-1111)

class HTauProcessor(ATLASStudent):

    def work(self):
                
        mc_tree = Tree(name = "_".join([self.fileset.name, "mc"]), model=TruthTau)
         
        reco_variables = (
            ("BDTJetScore", "F"),
            ("BDTEleScore", "F"),
            ("nPi0", "I"),
            ("pt", "F"),
            ("seedCalo_eta", "F"),
            ("seedCalo_phi", "F"),
            ("seedCalo_numTrack", "I"),
            ("charge", "I"),
        )

        truth_variables = (
            ("pt", "F"),
            ("m", "F"),
            ("eta", "F"),
            ("phi", "F"),
            ("vis_m", "F"),
            ("vis_eta", "F"),
            ("vis_phi", "F"),
            ("vis_Et", "F"),
            ("nProng", "I"),
            ("nPi0", "I"),
            ("charge", "I"),
        )

        common_variables = (
            ("matched", "I"),
            ("matched_dR", "F"),
            ("matched_collision", "B"),
        )

        # define branches for output ntuple
        variables = [
            ("Mvis_tau1_tau2","F"),
            ("numJets","I"),
            ("jetDeltaEta", "F"),
            ("numVertices","I"),
            ("MET","F"),
            ("MET_phi","F"),
            ("HT","F"),
            ("MMC_mass","F"),
            ("error", "B"),
            ("mu", "I"),
        ]

        if self.fileset.datatype == datasets.MC:
            variables.append(("selected", "B"))
        
        jet_variables = (
            ("E", "F"),
            ("pt", "F"),
            ("m", "F"),
            ("eta", "F"),
            ("phi", "F")
        )

        for v, t in reco_variables + common_variables:
            for tau in (1, 2):
                variables.append(("tau%i_%s" % (tau, v), t))
        
        # only create truth branches for MC
        if self.fileset.datatype != datasets.DATA:
            for v, t in truth_variables + common_variables:
                for tau in (1, 2):
                    variables.append(("trueTau%i_%s" % (tau, v), t))
        
        for v, t in jet_variables:
            for jet in (1, 2):
                variables.append(("jet%i_%s" % (jet, v), t))

        self.tree = TreeChain(self.fileset.treename, files = self.fileset.files)
        self.tree.init()
         
        copied_variables = []
        
        if self.fileset.datatype == datasets.MC:
            copied_variables = self.tree.glob("EF_*")+self.tree.glob("L1_*")+self.tree.glob("L2_*")
        
        buffer = TreeBuffer(variables)
        self.output.cd()
        self.D4PD = Tree(name = self.fileset.name)
        self.D4PD.set_branches_from_buffer(buffer)
        #self.D4PD.set_branches_from_buffer(self.tree.buffer, copied_variables, visible=False)
        
        # set the event filters
        # passthrough for MC for trigger acceptance studies
        if self.fileset.datatype != datasets.MC:
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
            self.tree.set_filters(self.event_filters)

        # define tree collections
        self.tree.define_collection(name="taus", prefix="tau_", size="tau_n")
        self.tree.define_collection(name="jets", prefix="jet_AntiKt4TopoEM_", size="jet_AntiKt4TopoEM_n")
        self.tree.define_collection(name="truetaus", prefix="trueTau_", size="trueTau_n")
        self.tree.define_collection(name="mc", prefix="mc_", size="mc_n", mixin=MCParticle)
        self.tree.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n")
        self.tree.define_collection(name="electrons", prefix="el_", size="el_n")
        self.tree.define_collection(name="vertices", prefix="vxp_", size="vxp_n")
         
        for event in self.tree:
             
            """
            Experimenting here....
            Need to get all MC tau final states to build ntuple for missing mass calculator pdfs
            Need to match jets to VBF jets
            Fill tree used for MMC
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
                        mc_tree.pt = decay.fourvect.Pt()
                        mc_tree.eta = decay.fourvect.Eta()
                        mc_tree.phi = decay.fourvect.Phi()
                        mc_tree.pt_vis = decay.fourvect_visible.Pt()
                        mc_tree.eta_vis = decay.fourvect_visible.Eta()
                        mc_tree.phi_vis = decay.fourvect_visible.Phi()
                        mc_tree.nu_pt = decay.fourvect_missing.Pt()
                        mc_tree.nu_eta = decay.fourvect_missing.Eta()
                        mc_tree.nu_phi = decay.fourvect_missing.Phi()
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
                if self.fileset.datatype == datasets.MC:
                    self.D4PD.selected = False
                    self.D4PD.Fill()
                continue
            if self.fileset.datatype == datasets.MC:
                self.D4PD.selected = True

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
            
            # HT
            sumET = event.MET_LocHadTopo_sumet + event.MET_MuonBoy_sumet - event.MET_RefMuon_Track_sumet
            self.D4PD.HT = sumET
            
            """
            Overlap removal
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
            MET
            """
            MET_LocHadTopo_etx = event.MET_LocHadTopo_etx_CentralReg + event.MET_LocHadTopo_etx_EndcapRegion + event.MET_LocHadTopo_etx_ForwardReg
            MET_LocHadTopo_ety = event.MET_LocHadTopo_ety_CentralReg + event.MET_LocHadTopo_ety_EndcapRegion + event.MET_LocHadTopo_ety_ForwardReg
            
            try:
                MET_MuonBoy_etx = event.MET_MuonBoy_etx_CentralReg + event.MET_MuonBoy_etx_EndcapRegion + event.MET_MuonBoy_etx_ForwardReg
                MET_MuonBoy_ety = event.MET_MuonBoy_ety_CentralReg + event.MET_MuonBoy_ety_EndcapRegion + event.MET_MuonBoy_ety_ForwardReg
            except AttributeError:
                MET_MuonBoy_etx = event.MET_MuonBoy_etx
                MET_MuonBoy_ety = event.MET_MuonBoy_ety

            try:
                MET_RefMuon_Track_etx = event.MET_RefMuon_Track_etx_CentralReg + event.MET_RefMuon_Track_etx_EndcapRegion + event.MET_RefMuon_Track_etx_ForwardReg
                MET_RefMuon_Track_ety = event.MET_RefMuon_Track_ety_CentralReg + event.MET_RefMuon_Track_ety_EndcapRegion + event.MET_RefMuon_Track_ety_ForwardReg
            except:
                MET_RefMuon_Track_etx = event.MET_RefMuon_Track_etx
                MET_RefMuon_Track_ety = event.MET_RefMuon_Track_ety
            
            METx = MET_LocHadTopo_etx + MET_MuonBoy_etx - MET_RefMuon_Track_etx
            METy = MET_LocHadTopo_ety + MET_MuonBoy_ety - MET_RefMuon_Track_ety

            MET = math.sqrt(METx**2 + METy**2)
            
            self.D4PD.MET = MET
            if MET > 0:
                phi = math.asin(METy / MET)
            else:
                phi = -1111.
            self.D4PD.MET_phi = phi
            
            """
            MMC and misc variables
            """
            self.D4PD.MMC_mass = missingMass(taus, jets, METx, METy, sumET, self.fileset.datatype)
            self.D4PD.Mvis_tau1_tau2 = utils.Mvis(taus[0].Et, taus[0].seedCalo_phi, taus[1].Et, taus[1].seedCalo_phi)
            self.D4PD.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or (vtx.type == 3 and vtx.nTracks >= 2)])
            self.D4PD.numJets = len(jets)
            if self.fileset.datatype == datasets.MC:
                self.D4PD.mu = event.lbn

            """
            Jet variables
            """
            # needs improvement (crack, JVF, etc...)
            forward_jets = [jet for jet in jets if jet.eta > 0]
            backward_jets = [jet for jet in jets if jet.eta < 0]

            best_forward_jet = None
            best_backward_jet = None
            
            if forward_jets:
                best_forward_jet = max(forward_jets, key=lambda jet: jet.E) 
            if backward_jets:
                best_backward_jet = max(backward_jets, key=lambda jet: jet.E)

            # forward jet is #1, backward jet is #2
            for i, jet in zip((1, 2), (best_forward_jet, best_backward_jet)):
                if jet:
                    for v, t in jet_variables:
                        setattr(self.D4PD, "jet%i_%s" % (i, v), getattr(jet, v))

            if best_forward_jet and best_backward_jet:
                self.D4PD.jetDeltaEta = best_forward_jet.eta - best_backward_jet.eta
            
            
            """
            Experimenting here....
            Need to match jets to VBF jets
            """ 
            if self.fileset.datatype == datasets.MC:
                if self.fileset.name.startswith("VBFH"):
                    VBF_partons = hepmc.get_VBF_partons(event)
            
            """
            Reco tau variables
            """
            for v, t in reco_variables:
                for i, tau in zip((1, 2), taus):
                    setattr(self.D4PD, "tau%i_%s" % (i, v), getattr(tau, v))
            
            """
            Truth-matching
            """
            if self.fileset.datatype == datasets.MC and len(event.truetaus) > 0:
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
                    self.D4PD.error = True
                matched_truth = []
                for i, tau in zip((1, 2), taus):
                    matching_truth_index = tau.trueTauAssoc_index
                    if matching_truth_index >= 0:
                        unmatched_reco.remove(i)
                        # check that this tau / true tau was not previously matched
                        if matching_truth_index+1 not in unmatched_truth or \
                           matching_truth_index in matched_truth:
                            print "ERROR: match collision!"
                            self.D4PD.tau1_matched_collision = True
                            self.D4PD.tau2_matched_collision = True
                            self.D4PD.trueTau1_matched_collision = True
                            self.D4PD.trueTau2_matched_collision = True
                            self.D4PD.error = True
                        else:
                            unmatched_truth.remove(matching_truth_index+1)
                            matched_truth.append(matching_truth_index)
                            setattr(self.D4PD, "tau%i_matched" % i, 1)
                            setattr(self.D4PD, "tau%i_matched_dR" % i, tau.trueTauAssoc_dr)
                            setattr(self.D4PD, "trueTau%i_matched" % i, 1)
                            setattr(self.D4PD, "trueTau%i_matched_dR" % i, event.truetaus[matching_truth_index].tauAssoc_dr)
                            for v, t in truth_variables:
                                setattr(self.D4PD, "trueTau%i_%s" % (i, v), getattr(event.truetaus[matching_truth_index], v))
                
                for i, j in zip(unmatched_reco, unmatched_truth):
                    for v, t in truth_variables:
                        setattr(self.D4PD, "trueTau%i_%s" % (i, v), getattr(event.truetaus[j-1], v))
             
            # fill output ntuple
            self.D4PD.Fill(reset=True)


def missingMass(taus, jets, METx, METy, sumET, datatype):
    """
    Missing mass calculation
    returns the most likely mass
    """    
    jetvec = ROOT.vector("TLorentzVector")()
    jetP4 = ROOT.TLorentzVector()
    jet_SumEt = 0.

    for jet in jets:
        jetP4.SetPtEtaPhiM(jet.pt/GeV, jet.eta, jet.phi, jet.m/GeV);
        jet_SumEt += jetP4.Et();
        jetvec.push_back(jetP4);

    MMC_SumEt = sumET/GeV
    MMC_SumEt -= taus[0].Et/GeV
    MMC_SumEt -= taus[1].Et/GeV
    MMC_SumEt -= jet_SumEt
    
    VisTau0 = ROOT.TLorentzVector()
    if taus[0].seedCalo_numTrack <= 1:
        tau0_decay_type = 10
        VisTau0.SetPtEtaPhiM(taus[0].pt/GeV, taus[0].seedCalo_eta, taus[0].seedCalo_phi, 0.8)
    else:
        tau0_decay_type = 30
        VisTau0.SetPtEtaPhiM(taus[0].pt/GeV, taus[0].seedCalo_eta, taus[0].seedCalo_phi, 1.2)
    
    VisTau1 = ROOT.TLorentzVector()
    if taus[1].seedCalo_numTrack <= 1:
        tau1_decay_type = 10
        VisTau1.SetPtEtaPhiM(taus[1].pt/GeV, taus[1].seedCalo_eta, taus[1].seedCalo_phi, 0.8)
    else:
        tau1_decay_type = 30
        VisTau1.SetPtEtaPhiM(taus[1].pt/GeV, taus[1].seedCalo_eta, taus[1].seedCalo_phi, 1.2)
    
    met_vec = ROOT.TVector2(METx/GeV, METy/GeV)
    
    mmc = ROOT.MissingMassCalculator()
    mmc.SetMetVec(met_vec)
    mmc.SetVisTauVec(0, VisTau0)
    mmc.SetVisTauVec(1, VisTau1)
    mmc.SetVisTauType(0, tau0_decay_type)
    mmc.SetVisTauType(1, tau1_decay_type)
    mmc.SetMetScanParamsUE(MMC_SumEt,  0., int(datatype == datasets.MC)) # data_type int : should be 1 for MC, 0 for data
    misMassTest = mmc.RunMissingMassCalculator()
    output_fitstatus = mmc.GetFitStatus() # MMC output: 1=found solution; 0= no slution
    MMC_mass = 0.
    if output_fitstatus == 1:
        MMC_mass = mmc.GetFittedMass(2) # use 2 instead of 1 to remove spikes in output
    return MMC_mass

import ROOT
import math

from rootpy.tree import Tree, TreeBuffer, TreeChain, TreeModel
from rootpy.math.physics.vector import LorentzVector, Vector3

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.batch import ATLASStudent

from higgstautau.mixins import MCParticle, TauFourMomentum
from higgstautau import tautools
from higgstautau.models import *

ROOT.gErrorIgnoreLevel = ROOT.kFatal

VERBOSE = False


class Event(TreeModel):

    match_collision = BoolCol(default=False)
    matched = BoolCol(False)
    MET_x = FloatCol()
    MET_y = FloatCol()
    MET_phi = FloatCol()
    MET = FloatCol()
    sumET = FloatCol()


class RecoTau(TreeModel):

    pt = FloatCol()
    phi = FloatCol()
    eta = FloatCol()
    fourvect = LorentzVector
    numTrack = IntCol()
    nPi0 = IntCol()

    privtx = Vector3
    privtx_chiSquared = FloatCol()
    privtx_numberDoF = FloatCol()
    privtx_jvf = FloatCol()

    secvtx = Vector3
    secvtx_chiSquared = FloatCol()
    secvtx_numberDoF = FloatCol()

    ipZ0SinThetaSigLeadTrk = FloatCol()
    ipSigLeadTrk = FloatCol()
    trFlightPathSig = FloatCol()


class MCTauProcessor(ATLASStudent):
    """
    This processor creates ntuples for missing mass calculator optimizations
    """

    def work(self):

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        tree = TreeChain(
                self.metadata.treename,
                files=self.files,
                events=self.events,
                cache=True,
                cache_size=10000000,
                learn_entries=30)

        self.output.cd()

        # this tree will contain info pertaining to true tau decays
        # for possible use in the optimization of a missing mass calculator
        out_tree = Tree(name="tau",
                model=(TrueTau_MCBlock.prefix('mctau1_') +
                       TrueTau_MCBlock.prefix('mctau2_') +
                       RecoTau.prefix('tau1_') +
                       RecoTau.prefix('tau2_') +
                       Event))

        tree.define_collection(
                name="mc", prefix="mc_", size="mc_n",
                mix=MCParticle)
        tree.define_collection(
                name="taus", prefix="tau_", size="tau_n",
                mix=TauFourMomentum)

        is_visible = lambda fourvect: (
                fourvect.Et() > 10 * GeV and abs(fourvect.Eta()) < 2.5)

        def closest_reco_tau(taus, mctau, dR=0.2):

            closest_tau = None
            closest_dR = 10000
            for tau in taus:
                dr = utils.dR(tau.eta, tau.phi, mctau.Eta(), mctau.Phi())
                if dr < dR and dr < closest_dR:
                    closest_tau = tau
                    closest_dR = dr
            return closest_tau

        for event in tree:

            """
            Need to get all MC tau final states to build ntuple for missing mass calculator pdfs
            """
            # Only accept taus from a Z or Higgs
            tau_decays = tautools.get_tau_decays(event, parent_pdgid=(23, 25))

            # There should be exactly two taus
            if len(tau_decays) != 2:
                #print "#MCTAUS != 2: %i" % len(tau_decays)
                """
                for decay in tautools.get_tau_decays(event, status=None):
                    print "status:"
                    print decay.init.status
                    print "parents:"
                    for parent in decay.init.iparents():
                        print parent.pdgId
                """
                continue
            decay1, decay2 = tau_decays

            if VERBOSE:
                print decay1
                print decay2

                print decay1.decay_length
                print decay2.decay_length

                print "%s -> %s" % (decay1.prod_vertex, decay1.decay_vertex)
                print "%s -> %s" % (decay2.prod_vertex, decay2.decay_vertex)
                print "===="

            nele1 = len(decay1.electrons)
            nmuon1 = len(decay1.muons)

            nele2 = len(decay2.electrons)
            nmuon2 = len(decay2.muons)

            out_tree.mctau1_visible = is_visible(decay1.fourvect_visible)
            out_tree.mctau1_electron = nele1 > 0
            out_tree.mctau1_muon = nmuon1 > 0
            out_tree.mctau1_hadronic = decay1.hadronic
            out_tree.mctau1_nprong = decay1.nprong
            out_tree.mctau1_npi0 = decay1.npi0
            out_tree.mctau1_nneutrals = decay1.nneutrals
            out_tree.mctau1_fourvect.set_from(decay1.fourvect)
            out_tree.mctau1_fourvect_vis.set_from(decay1.fourvect_visible)
            out_tree.mctau1_fourvect_miss.set_from(decay1.fourvect_missing)
            out_tree.mctau1_dR_tau_nu = decay1.dR_tau_nu
            out_tree.mctau1_dTheta3d_tau_nu = decay1.dTheta3d_tau_nu
            out_tree.mctau1_decay_length = decay1.decay_length
            out_tree.mctau1_prod_vertex.set_from(decay1.prod_vertex)
            out_tree.mctau1_decay_vertex.set_from(decay1.decay_vertex)

            out_tree.mctau2_visible = is_visible(decay2.fourvect_visible)
            out_tree.mctau2_electron = nele2 > 0
            out_tree.mctau2_muon = nmuon2 > 0
            out_tree.mctau2_hadronic = decay2.hadronic
            out_tree.mctau2_nprong = decay2.nprong
            out_tree.mctau2_npi0 = decay2.npi0
            out_tree.mctau2_nneutrals = decay2.nneutrals
            out_tree.mctau2_fourvect.set_from(decay2.fourvect)
            out_tree.mctau2_fourvect_vis.set_from(decay2.fourvect_visible)
            out_tree.mctau2_fourvect_miss.set_from(decay2.fourvect_missing)
            out_tree.mctau2_dR_tau_nu = decay2.dR_tau_nu
            out_tree.mctau2_dTheta3d_tau_nu = decay2.dTheta3d_tau_nu
            out_tree.mctau2_decay_length = decay2.decay_length
            out_tree.mctau2_prod_vertex.set_from(decay2.prod_vertex)
            out_tree.mctau2_decay_vertex.set_from(decay2.decay_vertex)

            # match to reco taus
            tau1 = closest_reco_tau(event.taus, decay1.fourvect_visible,
                                    dR=0.2)
            tau2 = closest_reco_tau(event.taus, decay2.fourvect_visible,
                                    dR=0.2)

            if tau1 is None or tau2 is None:
                # both taus not matched
                out_tree.matched = False
                out_tree.Fill()
                continue
            out_tree.matched = True

            if tau1 == tau2:
                # match collision
                out_tree.match_collision = True
                out_tree.Fill()
                continue
            out_tree.match_collision = False

            out_tree.tau1_pt = tau1.pt
            out_tree.tau1_phi = tau1.phi
            out_tree.tau1_eta = tau1.eta
            out_tree.tau1_fourvect.set_from(tau1.fourvect)
            out_tree.tau1_numTrack = tau1.numTrack
            out_tree.tau1_nPi0 = tau1.nPi0
            out_tree.tau1_privtx.set_from(Vector3(tau1.privtx_x,
                                                  tau1.privtx_y,
                                                  tau1.privtx_z))
            out_tree.tau1_privtx_chiSquared = tau1.privtx_chiSquared
            out_tree.tau1_privtx_numberDoF = tau1.privtx_numberDoF
            out_tree.tau1_privtx_jvf = tau1.privtx_jvf
            out_tree.tau1_secvtx.set_from(Vector3(tau1.secvtx_x,
                                                  tau1.secvtx_y,
                                                  tau1.secvtx_z))

            out_tree.tau2_pt = tau2.pt
            out_tree.tau2_phi = tau2.phi
            out_tree.tau2_eta = tau2.eta
            out_tree.tau2_fourvect.set_from(tau2.fourvect)
            out_tree.tau2_numTrack = tau2.numTrack
            out_tree.tau2_nPi0 = tau2.nPi0
            out_tree.tau2_privtx.set_from(Vector3(tau2.privtx_x,
                                                  tau2.privtx_y,
                                                  tau2.privtx_z))
            out_tree.tau2_privtx_chiSquared = tau2.privtx_chiSquared
            out_tree.tau2_privtx_numberDoF = tau2.privtx_numberDoF
            out_tree.tau2_secvtx.set_from(Vector3(tau2.secvtx_x,
                                                  tau2.secvtx_y,
                                                  tau2.secvtx_z))

            out_tree.tau1_ipZ0SinThetaSigLeadTrk = tau1.ipZ0SinThetaSigLeadTrk
            out_tree.tau1_ipSigLeadTrk = tau1.ipSigLeadTrk
            out_tree.tau1_trFlightPathSig = tau1.trFlightPathSig

            out_tree.tau2_ipZ0SinThetaSigLeadTrk = tau2.ipZ0SinThetaSigLeadTrk
            out_tree.tau2_ipSigLeadTrk = tau2.ipSigLeadTrk
            out_tree.tau2_trFlightPathSig = tau2.trFlightPathSig

            """
            tau_track_n
            tau_track_d0
            tau_track_z0
            tau_track_phi
            tau_track_theta
            tau_track_qoverp
            tau_track_pt
            tau_track_eta
            tau_track_atPV_d0
            tau_track_atPV_z0
            tau_track_atPV_phi
            tau_track_atPV_theta
            tau_track_atPV_qoverp
            tau_track_atPV_pt
            tau_track_atPV_eta
            tau_track_atTJVA_n
            tau_track_atTJVA_d0
            tau_track_atTJVA_z0
            tau_track_atTJVA_phi
            tau_track_atTJVA_theta
            tau_track_atTJVA_qoverp
            tau_track_atTJVA_pt
            tau_track_atTJVA_eta
            """

            if VERBOSE:
                print (tau1.privtx_x, tau1.privtx_y, tau1.privtx_z)
                print (tau1.secvtx_x, tau1.secvtx_y, tau1.secvtx_z)

            if self.metadata.year == 2011:
                out_tree.MET_x = event.MET_RefFinal_BDTMedium_etx
                out_tree.MET_y = event.MET_RefFinal_BDTMedium_ety
                out_tree.MET_phi = event.MET_RefFinal_BDTMedium_phi
                out_tree.MET = event.MET_RefFinal_BDTMedium_et
                out_tree.sumET = event.MET_RefFinal_BDTMedium_sumet
            else:
                out_tree.MET_x = event.MET_RefFinal_STVF_etx
                out_tree.MET_y = event.MET_RefFinal_STVF_ety
                out_tree.MET_phi = event.MET_RefFinal_STVF_phi
                out_tree.MET = event.MET_RefFinal_STVF_et
                out_tree.sumET = event.MET_RefFinal_STVF_sumet

            out_tree.Fill()

        self.output.cd()
        out_tree.FlushBaskets()
        out_tree.Write()

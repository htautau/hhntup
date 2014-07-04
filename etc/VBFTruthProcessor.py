import ROOT
import math
from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from rootpy.tree.filtering import *
from atlastools.filtering import GRLFilter
from higgstautau.filters import *
from higgstautau.filters import MCTriggersOLD as MCTriggers
from atlastools.batch import ATLASStudent
from rootpy.tree import Tree, TreeBuffer, TreeChain
from rootpy.tree.cutflow import Cutflow
from rootpy.vector import Vector2
from higgstautau.mixins import *
from higgstautau import hepmc
import higgstautau.tautools
from higgstautau.models import *
from higgstautau import eventshapes
from mmc import missingmass
from higgstautau import eventview

ROOT.gErrorIgnoreLevel = ROOT.kFatal


class VBFTruthProcessor(ATLASStudent):

    def work(self):

        D4PD_model = RecoTauBlock + RecoJetBlock + EventVariables + TrueTauBlock + PartonBlock

        # initialize the TreeChain of all input files (each containing one tree named self.metadata.treename)
        tree = TreeChain(self.metadata.treename,
                         files=self.files,
                         events=self.events,
                         cache=True,
                         cache_size=10000000,
                         learn_entries=30)

        # create output tree
        self.output.cd()
        D4PD = Tree(name=self.metadata.name, model=D4PD_model)

        copied_variables = ['actualIntPerXing',
                            'averageIntPerXing']
        copied_variables += mc_triggers

        D4PD.set_buffer(tree.buffer, variables=copied_variables, create_branches=True, visible=False)
        tree.always_read(copied_variables)

        # set the event filters
        # passthrough for MC for trigger acceptance studies
        self.event_filters = EventFilterList([
            GRLFilter(self.grl, passthrough = self.metadata.datatype != datasets.DATA),
            PriVertex(),
            LArError(),
            JetCleaningLoose(passthrough = self.metadata.datatype != datasets.DATA),
            #JetCleaningMedium(passthrough = self.metadata.datatype != datasets.DATA),
            LArHole(),
            #JetCrackVeto(),
            ElectronVeto(),
            MuonVeto(),
            TauElectronVeto(),
            TauMuonVeto(),
            TauAuthorTrack(),
        ])

        self.event_filters.insert(1, MCTriggers())
        tree.filters += self.event_filters

        cutflow = Cutflow()

        # define tree collections
        tree.define_collection(name="taus", prefix="tau_", size="tau_n", mix=TauFourMomentum)
        # jet_eta etc is AntiKt4LCTopo in tau-perf D3PDs
        tree.define_collection(name="jets", prefix="jet_", size="jet_n", mix=FourMomentum)
        tree.define_collection(name="truetaus", prefix="trueTau_", size="trueTau_n", mix=MCTauFourMomentum)
        tree.define_collection(name="mc", prefix="mc_", size="mc_n", mix=MCParticle)
        tree.define_collection(name="muons", prefix="mu_staco_", size="mu_staco_n")
        tree.define_collection(name="electrons", prefix="el_", size="el_n")
        tree.define_collection(name="vertices", prefix="vxp_", size="vxp_n")

        # define tree objects
        D4PD.define_object(name='tau1', prefix='tau1_')
        D4PD.define_object(name='tau2', prefix='tau2_')
        D4PD.define_object(name='jet1', prefix='jet1_')
        D4PD.define_object(name='jet2', prefix='jet2_')

        """
        tree.define_association(origin='taus', target='truetaus', prefix='trueTauAssoc_', link='index')
        tree.define_association(origin='truetaus', target='taus', prefix='tauAssoc_', link='index')
        """

        # entering the main event loop...
        for event in tree:

            D4PD.reset()
            cutflow.reset()

            # tau selection
            event.taus.select(lambda tau: tau.pt > 15*GeV)
            # Jet selection
            event.jets.select(lambda jet: jet.fourvect.P() > 25*GeV and abs(jet.emscale_eta) < 4.5)

            """
            Get VBF jets
            """
            # get partons (already sorted by eta in hepmc)
            parton1, parton2 = hepmc.get_VBF_partons(event)
            PartonBlock.set(D4PD, parton1, parton2)
            D4PD.dR_quarks = parton1.fourvect.DeltaR(parton2.fourvect)

            """
            Get true taus
            """
            event.truetaus.select(lambda tau: tau.vis_Et > 10 * GeV and abs(tau.vis_eta) < 2.5)
            if len(event.truetaus) > 2:
                print "ERROR: too many true taus: %i" % len(event.truetaus)
                D4PD.error = 1
                D4PD.Fill()
                continue
            elif len(event.truetaus) < 2:
                print "ERROR: too few true taus: %i" % len(event.truetaus)
                D4PD.error = 2
                D4PD.Fill()
                continue

            """
            fourvects = []
            colors = []
            radii = []
            for thing in event.jets:
                fourvects.append(thing.fourvect)
                colors.append('blue')
                radii.append(.4)
            for parton in (parton1, parton2):
                fourvects.append(parton.fourvect)
                colors.append('green')
                radii.append(.1)
            for thing in event.taus:
                fourvects.append(thing.fourvect)
                colors.append('red')
                radii.append(.2)
            for tau in event.truetaus:
                fourvects.append(tau.fourvect)
                colors.append('purple')
                radii.append(.1)

            eventview.draw(event, fourvects, colors=colors, radii=radii)
            """

            TrueTauBlock.set(D4PD, 1, event.truetaus[0])
            TrueTauBlock.set(D4PD, 2, event.truetaus[1])

            D4PD.dR_truetaus = event.truetaus[0].fourvect.DeltaR(event.truetaus[1].fourvect)

            D4PD.dR_quark_tau = min([
                                    parton1.fourvect.DeltaR(event.truetaus[0].fourvect),
                                    parton2.fourvect.DeltaR(event.truetaus[0].fourvect),
                                    parton1.fourvect.DeltaR(event.truetaus[1].fourvect),
                                    parton2.fourvect.DeltaR(event.truetaus[1].fourvect),
                                    ])

            taus = []
            if event.taus:
                for truetau in event.truetaus:
                    closest_tau = min(event.taus, key=lambda tau: utils.dR(tau.seedCalo_eta, tau.seedCalo_phi, truetau.vis_eta, truetau.vis_phi))
                    if utils.dR(closest_tau.seedCalo_eta, closest_tau.seedCalo_phi, truetau.vis_eta, truetau.vis_phi) < 0.2:
                        if closest_tau in taus:
                            # collision
                            D4PD.error = 3
                            break
                        taus.append(closest_tau)
            if len(taus) < 2:
                # collision
                D4PD.Fill()
                continue

            """
            # Overlap removal between taus and jets
            event.jets.select(lambda jet: not any([tau for tau in taus if
                              utils.dR(jet.emscale_eta, jet.emscale_phi, tau.seedCalo_eta, tau.seedCalo_phi) < .2]))
            """

            jets = []
            if event.jets:
                for quark in (parton1, parton2):
                    closest_jet = min(event.jets, key=lambda jet: utils.dR(jet.eta, jet.phi, quark.eta, quark.phi))
                    if utils.dR(closest_jet.eta, closest_jet.phi, quark.eta, quark.phi) < 0.4:
                        if closest_jet in jets:
                            # collision
                            D4PD.error = 4
                            break
                        jets.append(closest_jet)
            if len(jets) < 2:
                # collision
                D4PD.Fill()
                continue

            """
            Jet variables
            """
            RecoJetBlock.set(D4PD, jets[0], jets[1])

            """
            Reco tau variables
            This must come after the RecoJetBlock is filled since
            that sets the jet_beta for boosting the taus
            """
            RecoTauBlock.set(event, D4PD, taus[0], taus[1])
            D4PD.true_Mvis_tau1_tau2 = (D4PD.trueTau1_fourvect_vis + D4PD.trueTau2_fourvect_vis).M()

            """
            MET
            """
            METx = event.MET_RefFinal_etx
            METy = event.MET_RefFinal_ety
            MET_vect = Vector2(METx, METy)
            MET_3vect = Vector3(METx, METy, 0.)

            D4PD.MET = event.MET_RefFinal_et
            D4PD.MET_phi = event.MET_RefFinal_phi

            # HT TODO: Fix
            sumET = event.MET_RefFinal_sumet
            D4PD.HT = sumET

            D4PD.numVertices = len([vtx for vtx in event.vertices if (vtx.type == 1 and vtx.nTracks >= 4) or
                                    (vtx.type == 3 and vtx.nTracks >= 2)])
            # fill output ntuple
            # use reset=True to reset all variables to their defaults after the fill
            # to avoid any values from this event carrying over into the next
            D4PD.cutflow = cutflow.int()
            D4PD.Fill()

        self.output.cd()
        D4PD.FlushBaskets()
        D4PD.Write()

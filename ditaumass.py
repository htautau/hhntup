import ROOT
import math

from argparse import ArgumentParser

from rootpy.tree import Tree, TreeBuffer, TreeChain, TreeModel
from rootpy.math.physics.vector import LorentzVector, Vector3
from rootpy.hep import pdg

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.batch import ATLASStudent

from higgstautau.mixins import MCParticle, TauFourMomentum
from higgstautau import tautools
from higgstautau.models import *
from higgstautau.hadhad.objects import define_objects
from higgstautau.ditaumass.models import *
from higgstautau.utils import *

ROOT.gErrorIgnoreLevel = ROOT.kFatal


class ditaumass(ATLASStudent):

    def __init__(self, options, **kwargs):

        super(ditaumass, self).__init__(**kwargs)
        parser = ArgumentParser()
        parser.add_argument('--student-verbose',
            dest='verbose',
            action='store_true', default=False)
        parser.add_argument('--draw-decays',
            action='store_true', default=False)
        parser.add_argument('--higgs',
            action='store_true', default=False)
        self.args = parser.parse_args(options)

    def work(self):

        year = self.metadata.year
        verbose = self.args.verbose
        draw_decays = self.args.draw_decays
        args = self.args

        # initialize the TreeChain of all input files
        # only enable branches I need
        chain = TreeChain(
                self.metadata.treename,
                files=self.files,
                branches=[
                    'tau_*',
                    'mc_*',
                    'el_*',
                    'mu_staco_*',
                    'MET_RefFinal_BDTMedium_*',
                    'MET_RefFinal_STVF_*',
                    'EventNumber',
                    'RunNumber',
                    'averageIntPerXing',
                    ],
                events=self.events,
                read_branches_on_demand=True,
                cache=True,
                verbose=True)

        define_objects(chain, year)

        self.output.cd()

        # this tree will contain info pertaining to true tau decays
        # for possible use in the optimization of a missing mass calculator
        tree = Tree(name="ditaumass",
                model=(TrueTau.prefix('truetau1_') +
                       TrueTau.prefix('truetau2_') +
                       RecoTau.prefix('tau1_') +
                       RecoTau.prefix('tau2_') +
                       RecoElectron.prefix('ele1_') +
                       RecoElectron.prefix('ele2_') +
                       RecoMuon.prefix('muon1_') +
                       RecoMuon.prefix('muon2_') +
                       Event))

        tree.define_object(name='resonance', prefix='resonance_')
        tree.define_object(name='radiative', prefix='radiative_')

        truetaus = [
            tree.define_object(name='truetau1', prefix='truetau1_'),
            tree.define_object(name='truetau2', prefix='truetau2_')]

        taus = [
            tree.define_object(name='tau1', prefix='tau1_'),
            tree.define_object(name='tau2', prefix='tau2_')]

        electrons = [
            tree.define_object(name='ele1', prefix='ele1_'),
            tree.define_object(name='ele2', prefix='ele2_')]

        muons = [
            tree.define_object(name='muon1', prefix='muon1_'),
            tree.define_object(name='muon2', prefix='muon2_')]

        # get the Z or Higgs
        if args.higgs:
            resonance_pdgid = 25
        else:
            resonance_pdgid = 23

        if '7TeV' in self.metadata.name:
            collision_energy = 7
        else:
            collision_energy = 8

        for event_index, event in enumerate(chain):

            try:
                tree.reset_branch_values()

                # get the Z or Higgs
                resonance = tautools.get_particles(event, resonance_pdgid,
                        num_expected=1)

                if not resonance:
                    print "could not find resonance"
                    continue

                # get the resonance just before the decay
                resonance = resonance[0].last_self

                if draw_decays:
                    resonance.export_graphvis('resonance_%d.dot' %
                            event.EventNumber)

                FourVectModel.set(tree.resonance, resonance)

                # collect decay products (taus and photons)
                tau_decays = []
                mc_photons = []
                for child in resonance.iter_children():
                    if abs(child.pdgId) == pdg.tau_minus:
                        # ignore status 3 taus in 2012 (something strange in the
                        # MC record...)
                        if year == 2012:
                            if child.status == 3:
                                continue
                        tau_decays.append(tautools.TauDecay(child))
                    elif child.pdgId == pdg.gamma:
                        mc_photons.append(child)
                    else:
                        raise TypeError(
                                'unexpected particle after resonance:\n%s' %
                                child)

                # There should be exactly two taus
                if len(tau_decays) != 2:
                    print "found %i tau decays in MC record" % len(tau_decays)
                    for decay in tau_decays:
                        print decay
                    # skip this event
                    continue

                # check for incomplete tau decays
                invalid = False
                for decay in tau_decays:
                    if not decay.valid:
                        print "invalid tau decay:"
                        print decay
                        if draw_decays:
                            decay.init.export_graphvis(
                                    'decay_invalid_%d.dot' %
                                    event.EventNumber)
                        invalid = True
                if invalid:
                    # skip this event
                    continue

                radiative_fourvect = LorentzVector()
                for photon in mc_photons:
                    radiative_fourvect += photon.fourvect

                radiative_fourvect.fourvect = radiative_fourvect
                FourVectModel.set(tree.radiative, radiative_fourvect)
                tree.radiative_ngamma = len(mc_photons)
                tree.radiative_ngamma_5 = len([
                    ph for ph in mc_photons if ph.pt > 5])
                tree.radiative_ngamma_10 = len([
                    ph for ph in mc_photons if ph.pt > 10])
                tree.radiative_et_scalarsum = sum([
                    ph.pt for ph in mc_photons] + [0])

                all_matched = True
                matched_objects = []

                skip = False
                for i, (decay, truetau, tau, electron, muon) in enumerate(zip(
                        tau_decays, truetaus, taus, electrons, muons)):

                    if draw_decays:
                        decay.init.export_graphvis('decay%d_%d.dot' % (
                            i, event.EventNumber))

                    TrueTau.set(truetau, decay, verbose=verbose)

                    # match to reco taus, electrons and muons
                    if decay.hadronic:
                        recotau, dr = closest_reco_object(
                                event.taus, decay.fourvect_visible, dR=0.2)
                        if recotau is not None:
                            matched_objects.append(recotau)
                            recotau.matched = True
                            recotau.matched_dr = dr
                            RecoTau.set(tau, recotau, verbose=verbose)
                        else:
                            all_matched = False
                    elif decay.leptonic_electron:
                        recoele, dr = closest_reco_object(
                                event.electrons, decay.fourvect_visible, dR=0.2)
                        if recoele is not None:
                            matched_objects.append(recoele)
                            recoele.matched = True
                            recoele.matched_dr = dr
                            RecoElectron.set(electron, recoele)
                        else:
                            all_matched = False
                    elif decay.leptonic_muon:
                        recomuon, dr = closest_reco_object(
                                event.muons, decay.fourvect_visible, dR=0.2)
                        if recomuon is not None:
                            matched_objects.append(recomuon)
                            recomuon.matched = True
                            recomuon.matched_dr = dr
                            RecoMuon.set(muon, recomuon)
                        else:
                            all_matched = False
                    else:
                        print "unhandled invalid tau decay:"
                        print decay
                        if not draw_decays:
                            decay.init.export_graphvis('decay%d_%d.dot' % (
                                i, event.EventNumber))
                        # skip this event
                        skip = True
                        break
                if skip:
                    # skip this event
                    continue

                # did both decays match a reco object?
                tree.matched = all_matched

                # match collision: decays matched same reco object
                if all_matched:
                    tree.match_collision = (
                            matched_objects[0] == matched_objects[1])

                # MET
                tree.met_x = event.MET.etx
                tree.met_y = event.MET.ety
                tree.met_phi = event.MET.phi
                tree.met = event.MET.et
                tree.sum_et = event.MET.sumet

                # set extra event variables
                tree.channel = event.mc_channel_number
                tree.event = event.EventNumber
                tree.run = event.RunNumber
                tree.mu = event.averageIntPerXing
                tree.collision_energy = collision_energy

                tree.Fill()
            except:
                print "event index: %d" % event_index
                print "event number: %d" % event.EventNumber
                print "file: %s" % chain.file.GetName()
                raise

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

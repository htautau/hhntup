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

ROOT.gErrorIgnoreLevel = ROOT.kFatal

is_visible = lambda fourvect: (
                fourvect.Et() > 10 * GeV and abs(fourvect.Eta()) < 2.5)



class FourVectModel(TreeModel):

    pt = FloatCol()
    phi = FloatCol()
    eta = FloatCol()
    p = FloatCol()
    e = FloatCol()
    et = FloatCol()
    m = FloatCol()
    fourvect = LorentzVector

    @classmethod
    def set(cls, this, other):

        vect = other.fourvect
        this.pt = vect.Pt()
        this.p = vect.P()
        this.et = vect.Et()
        this.e = vect.E()
        this.m = vect.M()
        this.phi = vect.Phi()
        this.eta = vect.Eta()
        this.fourvect.set_from(vect)

    @classmethod
    def set_vis(cls, this, other):

        vect = other.fourvect_visible
        this.pt_vis = vect.Pt()
        this.p_vis = vect.P()
        this.et_vis = vect.Et()
        this.e_vis = vect.E()
        this.m_vis = vect.M()
        this.phi_vis = vect.Phi()
        this.eta_vis = vect.Eta()
        this.fourvect_vis.set_from(vect)

    @classmethod
    def set_miss(cls, this, other):

        vect = other.fourvect_missing
        this.pt_miss = vect.Pt()
        this.p_miss = vect.P()
        this.et_miss = vect.Et()
        this.e_miss = vect.E()
        this.m_miss = vect.M()
        this.phi_miss = vect.Phi()
        this.eta_miss = vect.Eta()
        this.fourvect_miss.set_from(vect)


class Event(FourVectModel.prefix('resonance_') +
            FourVectModel.prefix('radiative_')):

    match_collision = BoolCol(default=False)
    matched = BoolCol(False)

    met_x = FloatCol()
    met_y = FloatCol()
    met_phi = FloatCol()
    met = FloatCol()
    sum_et = FloatCol()

    radiative_ngamma = IntCol()


class TrueTau(FourVectModel +
        FourVectModel.suffix('_vis') +
        FourVectModel.suffix('_miss')):

    visible = BoolCol(default=False)
    charge = IntCol()

    hadronic = BoolCol(default=False)
    leptonic_electron = BoolCol(default=False)
    leptonic_muon = BoolCol(default=False)

    nprong = IntCol()
    npi_zero = IntCol()
    npi_ch = IntCol()
    nk_zero = IntCol()
    nk_ch = IntCol()
    nneutral = IntCol()
    nelectron = IntCol()
    nmuon = IntCol()
    nneutrino = IntCol()
    ngamma = IntCol()

    has_charged_rho = BoolCol()
    has_a1 = BoolCol()

    prod_vertex = Vector3
    decay_vertex = Vector3
    decay_length = FloatCol(default=-1111)

    dr_vistau_nu = FloatCol(default=-1111)
    dtheta3d_vistau_nu = FloatCol(default=-1111)

    @classmethod
    def set(cls, mctau, decay, verbose=False):

        mctau.visible = is_visible(decay.fourvect_visible)
        mctau.charge = decay.charge

        mctau.leptonic_electron = decay.leptonic_electron
        mctau.leptonic_muon = decay.leptonic_muon
        mctau.hadronic = decay.hadronic

        mctau.nprong = decay.nprong
        mctau.npi_zero = decay.npi0
        mctau.npi_ch = len(decay.charged_pions)
        mctau.nk_zero = len(decay.neutral_kaons)
        mctau.nk_ch = len(decay.charged_kaons)
        mctau.nneutral = decay.nneutrals
        mctau.nelectron = len(decay.electrons)
        mctau.nmuons = len(decay.muons)
        mctau.nneutrino = len(decay.neutrinos)
        mctau.ngamma = len(decay.photons)

        mctau.has_charged_rho = decay.has_charged_rho
        mctau.has_a1 = decay.has_a1

        FourVectModel.set(mctau, decay)
        FourVectModel.set_vis(mctau, decay)
        FourVectModel.set_miss(mctau, decay)

        mctau.dr_vistau_nu = decay.dr_vistau_nu
        mctau.dtheta3d_vistau_nu = decay.dtheta3d_vistau_nu

        mctau.decay_length = decay.decay_length
        mctau.prod_vertex.set_from(decay.prod_vertex)
        mctau.decay_vertex.set_from(decay.decay_vertex)

        # sanity check
        if mctau.hadronic:
            assert mctau.leptonic_electron == False
            assert mctau.leptonic_muon == False

        if verbose:
            print decay
            print decay.decay_length
            print "%s -> %s" % (decay.prod_vertex, decay.decay_vertex)
            print "===="


class RecoTau(FourVectModel):

    charge = IntCol()
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

    @classmethod
    def set(cls, tau, recotau, verbose=False):

        FourVectModel.set(tau, recotau)

        tau.charge = recotau.charge
        tau.numTrack = recotau.numTrack
        tau.nPi0 = recotau.nPi0
        tau.charge = recotau.charge

        tau.privtx.set_from(
           Vector3(recotau.privtx_x,
                   recotau.privtx_y,
                   recotau.privtx_z))
        tau.privtx_chiSquared = recotau.privtx_chiSquared
        tau.privtx_numberDoF = recotau.privtx_numberDoF
        tau.privtx_jvf = recotau.privtx_jvf

        tau.secvtx.set_from(
           Vector3(recotau.secvtx_x,
                   recotau.secvtx_y,
                   recotau.secvtx_z))
        tau.secvtx_chiSquared = recotau.secvtx_chiSquared
        tau.secvtx_numberDoF = recotau.secvtx_numberDoF

        tau.ipZ0SinThetaSigLeadTrk = recotau.ipZ0SinThetaSigLeadTrk
        tau.ipSigLeadTrk = recotau.ipSigLeadTrk
        tau.trFlightPathSig = recotau.trFlightPathSig

        if verbose:
            print
            print "reco tau:"
            print "privtx: ", (
                    recotau.privtx_x,
                    recotau.privtx_y,
                    recotau.privtx_z)
            print "secvtx: ", (
                    recotau.secvtx_x,
                    recotau.secvtx_y,
                    recotau.secvtx_z)
            print

        """ Additional variables to add later if needed
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


class RecoElectron(FourVectModel):

    charge = IntCol()

    @classmethod
    def set(cls, this, other):

        this.charge = other.charge
        FourVectModel.set(this, other)


class RecoMuon(FourVectModel):

    charge = IntCol()

    @classmethod
    def set(cls, this, other):

        this.charge = other.charge
        FourVectModel.set(this, other)


def closest_reco_object(objects, thing, dR=0.2):

    closest_object = None
    closest_dR = 10000
    for other in objects:
        dr = utils.dR(other.eta, other.phi, thing.Eta(), thing.Phi())
        if dr < dR and dr < closest_dR:
            closest_object = other
            closest_dR = dr
    return closest_object


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
                    ],
                events=self.events,
                cache=True,
                cache_size=10000000,
                learn_entries=30,
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

        for event_index, event in enumerate(chain):

            try:
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
                    # skip event
                    continue

                # check for incomplete tau decays
                incomplete = False
                for decay in tau_decays:
                    if not decay.complete:
                        print "found incomplete tau decay:\n%s" % decay
                        incomplete = True
                if incomplete:
                    # skip event
                    continue

                radiative_fourvect = LorentzVector()
                for photon in mc_photons:
                    radiative_fourvect += photon.fourvect

                radiative_fourvect.fourvect = radiative_fourvect
                FourVectModel.set(tree.radiative, radiative_fourvect)
                tree.radiative_ngamma = len(mc_photons)

                matched = True
                matched_objects = []

                for i, (decay, truetau, tau, electron, muon) in enumerate(zip(
                        tau_decays, truetaus, taus, electrons, muons)):

                    if draw_decays:
                        decay.init.export_graphvis('decay%d_%d.dot' % (
                            i, event.EventNumber))

                    TrueTau.set(truetau, decay, verbose=verbose)

                    # match to reco taus, electrons and muons
                    if decay.hadronic:
                        recotau = closest_reco_object(
                                event.taus, decay.fourvect_visible, dR=0.2)
                        matched_objects.append(recotau)
                        if recotau is not None:
                            RecoTau.set(tau, recotau, verbose=verbose)
                        else:
                            matched = False
                    elif decay.leptonic_electron:
                        recoele = closest_reco_object(
                                event.electrons, decay.fourvect_visible, dR=0.2)
                        matched_objects.append(recoele)
                        if recoele is not None:
                            RecoElectron.set(electron, recoele)
                        else:
                            matched = False
                    elif decay.leptonic_muon:
                        recomuon = closest_reco_object(
                                event.muons, decay.fourvect_visible, dR=0.2)
                        matched_objects.append(recomuon)
                        if recomuon is not None:
                            RecoMuon.set(muon, recomuon)
                        else:
                            matched = False
                    else:
                        print decay
                        raise TypeError("Invalid tau decay")

                # did both decays match a reco object?
                tree.matched = matched

                # match collision: decays matched same reco object
                tree.match_collision = matched_objects[0] == matched_objects[1]

                # MET
                tree.met_x = event.MET.etx
                tree.met_y = event.MET.ety
                tree.met_phi = event.MET.phi
                tree.met = event.MET.et
                tree.sum_et = event.MET.sumet

                tree.Fill(reset=True)
            except:
                print "event index: %d" % event_index
                print "event number: %d" % event.EventNumber
                print "file: %s" % chain.file.GetName()
                raise

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

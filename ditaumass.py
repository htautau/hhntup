import ROOT
import math

from argparse import ArgumentParser

from rootpy.tree import Tree, TreeBuffer, TreeChain, TreeModel
from rootpy.math.physics.vector import LorentzVector, Vector3

from atlastools import datasets
from atlastools import utils
from atlastools.units import *
from atlastools.batch import ATLASStudent

from higgstautau.mixins import MCParticle, TauFourMomentum
from higgstautau import tautools
from higgstautau.models import *
from higgstautau.hadhad.objects import define_objects

ROOT.gErrorIgnoreLevel = ROOT.kFatal

VERBOSE = True
EXPORT_GRAPHVIZ = False

is_visible = lambda fourvect: (
                fourvect.Et() > 10 * GeV and abs(fourvect.Eta()) < 2.5)


class Event(TreeModel):

    match_collision = BoolCol(default=False)
    matched = BoolCol(False)
    MET_x = FloatCol()
    MET_y = FloatCol()
    MET_phi = FloatCol()
    MET = FloatCol()
    sumET = FloatCol()


class FourVectModel(TreeModel):

    charge = IntCol()
    pt = FloatCol()
    phi = FloatCol()
    eta = FloatCol()
    p = FloatCol()
    fourvect = LorentzVector

    @classmethod
    def set(cls, this, other):

        this.charge = other.charge
        vect = other.fourvect
        this.pt = vect.Pt()
        this.phi = vect.Phi()
        this.eta = vect.Eta()
        this.p = vect.P()
        this.fourvect.set_from(vect)

    @classmethod
    def set_vis(cls, this, other):

        vect = other.fourvect_visible
        this.pt_vis = vect.Pt()
        this.phi_vis = vect.Phi()
        this.eta_vis = vect.Eta()
        this.p_vis = vect.P()
        this.fourvect_vis.set_from(vect)

    @classmethod
    def set_miss(cls, this, other):

        vect = other.fourvect_missing
        this.pt_miss = vect.Pt()
        this.phi_miss = vect.Phi()
        this.eta_miss = vect.Eta()
        this.p_miss = vect.P()
        this.fourvect_miss.set_from(vect)


class TrueTau(FourVectModel +
        FourVectModel.suffix('_vis') +
        FourVectModel.suffix('_miss')):

    visible = BoolCol(default=False)
    hadronic = BoolCol(default=False)
    electron = BoolCol(default=False)
    muon = BoolCol(default=False)
    nprong = IntCol(default=-1111)
    npi0 = IntCol(default=-1111)
    nneutrals = IntCol(default=-1111)

    prod_vertex = Vector3
    decay_vertex = Vector3
    decay_length = FloatCol(default=-1111)

    dR_vistau_nu = FloatCol(default=-1111)
    dTheta3d_vistau_nu = FloatCol(default=-1111)

    @classmethod
    def set(cls, mctau, decay):

        mctau.visible = is_visible(decay.fourvect_visible)
        mctau.electron = decay.electron
        mctau.muon = decay.muon
        mctau.hadronic = decay.hadronic
        mctau.nprong = decay.nprong
        mctau.npi0 = decay.npi0
        mctau.nneutrals = decay.nneutrals
        mctau.charge = decay.charge

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
            assert mctau.electron == False
            assert mctau.muon == False

        if VERBOSE:
            print decay
            print decay.decay_length
            print "%s -> %s" % (decay.prod_vertex, decay.decay_vertex)
            print "===="


class RecoTau(FourVectModel):

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
    def set(cls, tau, recotau):

        FourVectModel.set(tau, recotau)

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

        tau.ipZ0SinThetaSigLeadTrk = recotau.ipZ0SinThetaSigLeadTrk
        tau.ipSigLeadTrk = recotau.ipSigLeadTrk
        tau.trFlightPathSig = recotau.trFlightPathSig

        if VERBOSE:
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

    pass

class RecoMuon(FourVectModel):

    pass


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
        parser.add_argument('--export-graphviz',
            action='store_true', default=False)
        self.args = parser.parse_args(options)
        global VERBOSE
        VERBOSE = self.args.verbose
        global EXPORT_GRAPHVIZ
        EXPORT_GRAPHVIZ = self.args.export_graphviz

    def work(self):

        year = self.metadata.year

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

        for event in chain:

            # Only accept taus from a Z or Higgs
            tau_decays = tautools.get_tau_decays(event,
                    parent_pdgid=(23, 25),
                    num_expected=2)

            # There should be exactly two taus
            if len(tau_decays) != 2:
                print "#MCTAUS != 2: %i" % len(tau_decays)
                for decay in tautools.get_tau_decays(event, status=None):
                    print "status:"
                    print decay.init.status
                    print "parents:"
                    for parent in decay.init.iter_parents():
                        print parent.pdgId
                # skip this event
                continue

            matched = True
            matched_objects = []

            for i, (decay, truetau, tau, electron, muon) in enumerate(zip(
                    tau_decays, truetaus, taus, electrons, muons)):

                if EXPORT_GRAPHVIZ:
                    decay.init.export_graphvis('decay%d_%d.dot' % (
                        i, event.EventNumber))

                TrueTau.set(truetau, decay)

                # match to reco taus, electrons and muons
                if decay.hadronic:
                    recotau = closest_reco_object(
                            event.taus, decay.fourvect_visible, dR=0.2)
                    matched_objects.append(recotau)
                    if recotau is not None:
                        RecoTau.set(tau, recotau)
                    else:
                        matched = False
                elif decay.electron:
                    recoele = closest_reco_object(
                            event.electrons, decay.fourvect_visible, dR=0.2)
                    matched_objects.append(recoele)
                    if recoele is not None:
                        RecoElectron.set(electron, recoele)
                    else:
                        matched = False
                elif decay.muon:
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
            tree.MET_x = event.MET.etx
            tree.MET_y = event.MET.ety
            tree.MET_phi = event.MET.phi
            tree.MET = event.MET.et
            tree.sumET = event.MET.sumet

            tree.Fill()

        self.output.cd()
        tree.FlushBaskets()
        tree.Write()

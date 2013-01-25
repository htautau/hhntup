from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzVector
from rootpy.types import *


class MatchedObject(TreeModel):

    matched = BoolCol()
    match_dr = FloatCol(default=1111)


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

    collision_energy = FloatCol()

    event = IntCol()
    run = IntCol()
    channel = IntCol()
    mu = FloatCol()

    match_collision = BoolCol(default=False)
    matched = BoolCol(False)

    met_x = FloatCol()
    met_y = FloatCol()
    met_phi = FloatCol()
    met = FloatCol()
    sum_et = FloatCol()

    radiative_ngamma = IntCol()
    radiative_ngamma_5 = IntCol()
    radiative_ngamma_10 = IntCol()
    radiative_et_scalarsum = FloatCol()


class TrueTau(FourVectModel +
        FourVectModel.suffix('_vis') +
        FourVectModel.suffix('_miss'),
        MatchedObject):

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


class RecoTau(FourVectModel, MatchedObject):

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

    dr_nu = FloatCol(default=-1111)
    dtheta3d_nu = FloatCol(default=-1111)

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


class RecoElectron(FourVectModel, MatchedObject):

    charge = IntCol()

    @classmethod
    def set(cls, this, other):

        this.charge = other.charge
        FourVectModel.set(this, other)


class RecoMuon(FourVectModel, MatchedObject):

    charge = IntCol()

    @classmethod
    def set(cls, this, other):

        this.charge = other.charge
        FourVectModel.set(this, other)

import ROOT
from rootpy.tree import TreeModel
from rootpy.math.physics.vector import LorentzVector, Vector3
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


class DecayVertex(TreeModel):

    privtx_x = FloatCol()
    privtx_y = FloatCol()
    privtx_z = FloatCol()

    secvtx_x = FloatCol()
    secvtx_y = FloatCol()
    secvtx_z = FloatCol()

    privtx = Vector3
    secvtx = Vector3

    decay_theta = FloatCol()
    decay_phi = FloatCol()
    decay_length = FloatCol(default=-1111)

    @classmethod
    def set(cls, outobj, inobj):

        outobj.privtx_x = inobj.privtx_x
        outobj.privtx_y = inobj.privtx_y
        outobj.privtx_z = inobj.privtx_z

        outobj.secvtx_x = inobj.secvtx_x
        outobj.secvtx_y = inobj.secvtx_y
        outobj.secvtx_z = inobj.secvtx_z

        outobj.privtx = inobj.privtx
        outobj.secvtx = inobj.secvtx

        outobj.decay_theta = inobj.decay_vect.Theta()
        outobj.decay_phi = inobj.decay_vect.Phi()
        outobj.decay_length = inobj.decay_vect.Mag()


class RecoDecayVertex(DecayVertex):

    # covariance matrix
    privtx_xx = FloatCol()
    privtx_yy = FloatCol()
    privtx_zz = FloatCol()
    privtx_xy = FloatCol()
    privtx_yz = FloatCol()
    privtx_zx = FloatCol()

    secvtx_xx = FloatCol()
    secvtx_yy = FloatCol()
    secvtx_zz = FloatCol()
    secvtx_xy = FloatCol()
    secvtx_yz = FloatCol()
    secvtx_zx = FloatCol()

    privtx_chiSquared = FloatCol()
    privtx_numberDoF = FloatCol()
    privtx_jvf = FloatCol()
    # ROOT.TMath.Prob(privtx_chiSquared, privtx_numberDoF)
    privtx_prob = FloatCol()

    secvtx_chiSquared = FloatCol()
    secvtx_numberDoF = FloatCol()
    # ROOT.TMath.Prob(secvtx_chiSquared, secvtx_numberDoF)
    secvtx_prob = FloatCol()

    decay_length_sigma = FloatCol()
    # L / L_sigma
    decay_length_significance = FloatCol()

    @classmethod
    def set(cls, outobj, inobj):

        DecayVertex.set(outobj, inobj)

        outobj.privtx_xx = inobj.privtx_xx
        outobj.privtx_yy = inobj.privtx_yy
        outobj.privtx_zz = inobj.privtx_zz
        outobj.privtx_xy = inobj.privtx_xy
        outobj.privtx_yz = inobj.privtx_yz
        outobj.privtx_zx = inobj.privtx_zx

        outobj.secvtx_xx = inobj.secvtx_xx
        outobj.secvtx_yy = inobj.secvtx_yy
        outobj.secvtx_zz = inobj.secvtx_zz
        outobj.secvtx_xy = inobj.secvtx_xy
        outobj.secvtx_yz = inobj.secvtx_yz
        outobj.secvtx_zx = inobj.secvtx_zx

        outobj.privtx_chiSquared = inobj.privtx_chiSquared
        outobj.privtx_numberDoF = inobj.privtx_numberDoF
        outobj.privtx_jvf = inobj.privtx_jvf
        outobj.privtx_prob = ROOT.TMath.Prob(
                inobj.privtx_chiSquared,
                int(inobj.privtx_numberDoF))

        outobj.secvtx_chiSquared = inobj.secvtx_chiSquared
        outobj.secvtx_numberDoF = inobj.secvtx_numberDoF
        outobj.secvtx_prob = ROOT.TMath.Prob(
                inobj.secvtx_chiSquared,
                int(inobj.secvtx_numberDoF))

        x = inobj.decay_vect.X()
        y = inobj.decay_vect.Y()
        z = inobj.decay_vect.Z()

        dx2 = inobj.privtx_xx
        dy2 = inobj.privtx_yy
        dz2 = inobj.privtx_zz

        decay_length_sigma =
        decay_length_significance = inobj.decay_length / decay_length_sigma

        outobj.decay_length_sigma = decay_length_sigma
        outobj.decay_length_significance = decay_length_significance


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
        MatchedObject,
        DecayVertex):

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
        mctau.privtx.set_from(decay.prod_vertex)
        mctau.secvtx.set_from(decay.decay_vertex)

        # sanity check
        if mctau.hadronic:
            assert mctau.leptonic_electron == False
            assert mctau.leptonic_muon == False

        if verbose:
            print decay
            print decay.decay_length
            print "%s -> %s" % (decay.prod_vertex, decay.decay_vertex)
            print "===="


class RecoTau(FourVectModel, MatchedObject, RecoDecayVertex):

    charge = IntCol()
    numTrack = IntCol()
    nPi0 = IntCol()

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

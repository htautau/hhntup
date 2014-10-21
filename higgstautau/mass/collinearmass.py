import math
import ROOT
from rootpy.vector import LorentzVector


def is_MET_bisecting(dphi_taus, dphi_tau1_MET, dphi_tau2_MET):
    """
    check whether the MET is bisecting the
    two taus in the transverse plane

    dphi_taus = between the 2 taus
    dphi_tau1_MET = between tau 1 and MET
    dphi_tau2_MET = between tau 2 and MET
    """
    return (((max(dphi_tau1_MET, dphi_tau2_MET) <= dphi_taus) and
            (dphi_tau1_MET + dphi_tau2_MET <= math.pi)))


def mass(tau1, tau2, METpx, METpy):
    """
    Calculate and return the collinear mass and momentum fractions
    of tau1 and tau2

    TODO: set visible mass of taus. 1.2 GeV for 3p and 0.8 GeV for 1p
    """
    recTau1 = LorentzVector()
    recTau2 = LorentzVector()

    # tau 4-vector; synchronize for MMC calculation
    if tau1.obj.nTracks() < 3:
        recTau1.SetPtEtaPhiM(tau1.obj.pt(), tau1.obj.eta(), tau1.obj.phi(), 800.) # MeV
    else:
        recTau1.SetPtEtaPhiM(tau1.obj.pt(), tau1.obj.eta(), tau1.obj.phi(), 1200.) # MeV

    if tau2.obj.nTracks() < 3:
        recTau2.SetPtEtaPhiM(tau2.obj.pt(), tau2.obj.eta(), tau2.obj.phi(), 800.) # MeV
    else:
        recTau2.SetPtEtaPhiM(tau2.obj.pt(), tau2.obj.eta(), tau2.obj.phi(), 1200.) # MeV

    K = ROOT.TMatrixD(2, 2)
    K[0][0] = recTau1.Px(); K[0][1] = recTau2.Px()
    K[1][0] = recTau1.Py(); K[1][1] = recTau2.Py()

    if K.Determinant() == 0:
        return -1., -1111., -1111.

    M = ROOT.TMatrixD(2, 1)
    M[0][0] = METpx
    M[1][0] = METpy

    Kinv = K.Invert()

    X = Kinv * M

    X1 = X(0, 0)
    X2 = X(1, 0)

    x1 = 1./(1. + X1)
    x2 = 1./(1. + X2)

    p1 = recTau1 * (1. / x1)
    p2 = recTau2 * (1. / x2)
    m_col = (p1 + p2).M()
    m_vis = (recTau1 + recTau2).M()

    return m_vis, m_col, x1, x2

import math
import ROOT
from rootpy.math.physics.vector import LorentzVector


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
    if tau1.numTrack < 3:
        recTau1.SetPtEtaPhiM(tau1.pt, tau1.eta, tau1.phi, 800.) # MeV
    else:
        recTau1.SetPtEtaPhiM(tau1.pt, tau1.eta, tau1.phi, 1200.) # MeV

    if tau2.numTrack<3:
        recTau2.SetPtEtaPhiM(tau2.pt, tau2.eta, tau2.phi, 800.) # MeV
    else:
        recTau2.SetPtEtaPhiM(tau2.pt, tau2.eta, tau2.phi, 1200.) # MeV

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


def mass_soshi(tau1, tau2, METpx, METpy):
    """
    Collinear Approximation
    Code by Soshi Tsuno
    """
    recTau1 = LorentzVector()
    recTau2 = LorentzVector()

    # tau 4-vector; synchronize for MMC calculation
    if tau1.numTrack < 3:
        recTau1.SetPtEtaPhiM(tau1.pt, tau1.eta, tau1.phi, 800.) # MeV
    else:
        recTau1.SetPtEtaPhiM(tau1.pt, tau1.eta, tau1.phi, 1200.) # MeV

    if tau2.numTrack<3:
        recTau2.SetPtEtaPhiM(tau2.pt, tau2.eta, tau2.phi, 800.) # MeV
    else:
        recTau2.SetPtEtaPhiM(tau2.pt, tau2.eta, tau2.phi, 1200.) # MeV

    #recMet.SetPxPyPzE(RecalcMET_etx/1000.0,RecalcMET_ety/1000.0,0.0,RecalcMET_et/1000.0)

    recX1 = -1111.
    recX2 = -1111.

    denomRec1 = (recTau2.Py() * METpx - recTau2.Px() * METpy + recTau2.Py() * recTau1.Px() - recTau2.Px() * recTau1.Py())
    denomRec2 = (recTau1.Py() * METpx - recTau1.Px() * METpy + recTau1.Py() * recTau2.Px() - recTau1.Px() * recTau2.Py())

    if denomRec1 != 0.:
        recX1 = (recTau2.Py() * recTau1.Px() - recTau2.Px() * recTau1.Py()) / denomRec1
    if denomRec2 != 0.:
        recX2 = (recTau1.Py() * recTau2.Px() - recTau1.Px() * recTau2.Py()) / denomRec2

    recTauCol1 = recTau1.Clone()
    recTauCol2 = recTau2.Clone()
    if recX1 > -1111.:
        recTauCol1 *= 1. / recX1
    if recX2 > -1111.:
        recTauCol2 *= 1. / recX2

    recVis = recTau1 + recTau2
    recCol = recTauCol1 + recTauCol2

    return recVis.M(), recCol.M(), recX1, recX2

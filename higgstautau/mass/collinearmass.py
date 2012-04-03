import math
import ROOT

def is_MET_bisecting(dphi, dphi1, dphi2):
    """
    check whether the MET is bisecting the
    two taus in the transverse plane

    dphi = between the 2 taus
    dphi1 = between tau 1 and MET
    dphi2 = beteen tau 2 and MET
    """
    return ((max(dphi1, dphi2) <= dphi) and (dphi1 + dphi2 <= math.pi))


def mass(tau1, tau2, METpx, METpy):
    """
    Calculate and return the collinear mass and momentum fractions
    of tau1 and tau2
    """
    K = ROOT.TMatrixD(2, 2)
    K[0][0] = tau1.fourvect.Px(); K[0][1] = tau2.fourvect.Px()
    K[1][0] = tau1.fourvect.Py(); K[1][1] = tau2.fourvect.Py()

    if K.Determinant() == 0:
        return -1, -1111, -1111

    M = ROOT.TMatrixD(2, 1)
    M[0][0] = METpx
    M[1][0] = METpy

    Kinv = K.Invert()

    X = Kinv * M

    X1 = X(0, 0)
    X2 = X(1, 0)

    x1 = 1./(1. + X1)
    x2 = 1./(1. + X2)

    p1 = tau1.fourvect * (1. / x1)
    p2 = tau2.fourvect * (1. / x2)
    m = (p1 + p2).M()

    return m, x1, x2

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

    TODO: set visible mass of taus. 1.2 GeV for 3p and 0.8 GeV for 1p
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


"""
//=========================
// Collinear Approximation
//=========================
double recX1 = -100.0;
double recX2 = -100.0;
TLorentzVector recTau1,recTau2,recMet;

// tau 4-vector try to synchronize for MMC calculation
if (tau1_numTrack<3) { recTau1.SetPtEtaPhiM(tau1_pt/1000.0,tau1_eta,tau1_phi,0.8); }
else                 { recTau1.SetPtEtaPhiM(tau1_pt/1000.0,tau1_eta,tau1_phi,1.2); }

if (tau2_numTrack<3) { recTau2.SetPtEtaPhiM(tau2_pt/1000.0,tau2_eta,tau2_phi,0.8); }
else                 { recTau2.SetPtEtaPhiM(tau2_pt/1000.0,tau2_eta,tau2_phi,1.2); }

recMet.SetPxPyPzE(RecalcMET_etx/1000.0,RecalcMET_ety/1000.0,0.0,RecalcMET_et/1000.0);

double denomRec1 = (recTau2.Py()*recMet.Px()-recTau2.Px()*recMet.Py()+recTau2.Py()*recTau1.Px()-recTau2.Px()*recTau1.Py());
double denomRec2 = (recTau1.Py()*recMet.Px()-recTau1.Px()*recMet.Py()+recTau1.Py()*recTau2.Px()-recTau1.Px()*recTau2.Py());
if (denomRec1 != 0.0) { recX1 = (recTau2.Py()*recTau1.Px()-recTau2.Px()*recTau1.Py())/denomRec1; }
if (denomRec2 != 0.0) { recX2 = (recTau1.Py()*recTau2.Px()-recTau1.Px()*recTau2.Py())/denomRec2; }

TLorentzVector recTauCol1 = recTau1;
TLorentzVector recTauCol2 = recTau2;
if (recX1 > -100.0) { recTauCol1 *= 1.0/recX1; }
if (recX2 > -100.0) { recTauCol2 *= 1.0/recX2; }

TLorentzVector recVis = recTau1 + recTau2;
TLorentzVector recCol = recTauCol1 + recTauCol2;
"""

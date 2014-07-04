import ROOT
import os

ROOT.gSystem.Load(os.path.join(os.path.dirname(__file__), "libditaumass.so"))

from ROOT import DiTauMass
HAD1P = DiTauMass.Modes.HAD1P
HAD3P = DiTauMass.Modes.HAD3P
LEP = DiTauMass.Modes.LEP

MET_NSIGMA = 3.
# HWHM of Cauchy distribution (as fraction of tau mass)
TAU_MASS_HWHM = 0.05
USE_TNC = True
VERBOSE = 0

calculator = DiTauMass.Calculator(MET_NSIGMA, TAU_MASS_HWHM, VERBOSE, USE_TNC)


def mass(*args):

    return calculator.mass(*args)


def mass_scan(*args):

    return calculator.mass_scan(*args)

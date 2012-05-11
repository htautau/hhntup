#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Compile and load the egamma EnergyRescaler correction tool
print "    -----Loading the egamma EnergyRescaler correction tool (and compiling if necessary)..."
compileC(filedir + "/EnergyRescaler.cxx")
from ROOT import AnalysisFramework
egammaER = AnalysisFramework.External.EnergyRescaler()
egammaER.useDefaultCalibConstants("2011")

# Compile and load the egamma CaloIsoCorrection tool
print "    -----Loading the egamma CaloIsoCorrection tool (and compiling if necessary)..."
compileC(filedir + "/CaloIsoCorrection.cxx")
egammaCIC = AnalysisFramework.External.CaloIsoCorrection()

# Compile and load the egamma SF correction tool
print "    -----Loading the egamma SF correction tool (and compiling if necessary)..."
compileC(filedir + "/egammaSFclass.cxx")
egammaSF = AnalysisFramework.External.egammaSFclass()


#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ROOT import vector
from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))
__cached__ = False

# Compile and load the AnalysisMuonEfficiencyScaleFactors correction tool
print "    -----Loading the AnalysisMuonEfficiencyScaleFactors correction tool..."
compileC(filedir + "/EtaPhiBinning.cxx")
compileC(filedir + "/MuonEfficiencyScaleFactor.h")
compileC(filedir + "/AnalysisMuonEfficiencyScaleFactors.cxx")
from ROOT import AnalysisFramework
sharedir = filedir + "/SF_share"
def getMuonSF(pileup):
    global __cached__
    if not __cached__:
        int_lum = pileup.getIntegratedLumiVector()
        __cached__ = AnalysisFramework.External.AnalysisMuonEfficiencyScaleFactors("STACO_CB", int_lum, "MeV", sharedir)
    return __cached__

# Compile and load the muon SmearingClass tool
print "    -----Loading the muon SmearingClass tool (and compiling if necessary)..."
compileC(filedir + "/SmearingClass.cxx")
muonSmear = AnalysisFramework.External.SmearingClass("Data11", "staco", "pT", "Rel17", filedir + "/Smear_share/")

# Compile and load the muon trigger SF tool
print "    -----Loading the muon trigger SF tool (and compiling if necessary)..."
compileC(filedir + "/LeptonTriggerSF.cxx")
muonTriggerSF = AnalysisFramework.External.LeptonTriggerSF(filedir+"/Trigger_share/")

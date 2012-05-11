#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Compile and load the Moriond MMC tool
print "    -----Loading Moriond's MissingMassCalculator tool (and compiling if necessary)..."
compileC(filedir+"/Moriond/JERProviderMMC_Moriond.cxx")
compileC(filedir+"/Moriond/MissingMassCalculator_Moriond.cxx")

from ROOT import AnalysisFramework
mmcTool_moriond = AnalysisFramework.External.MissingMassCalculator(filedir+"/Moriond/JERProviderPlots.root")
mmcTool_moriond.SetAlgorithmVersion(1)
mmcTool_moriond.SetNiterFit1(20)
mmcTool_moriond.SetNiterFit2(30)
mmcTool_moriond.SetNiterFit3(10)
mmcTool_moriond.SetNsigmaMETscan(3.0)
mmcTool_moriond.SetApplyMassScale(0)

## Compile and load the MMC tool
#print "    -----Loading the MissingMassCalculator tool (and compiling if necessary)..."
#compileC(filedir+"/JERProviderMMC.cxx")
#compileC(filedir+"/MissingMassCalculator.cxx")

#from ROOT import AnalysisFramework
#mmcTool = AnalysisFramework.External.MissingMassCalculator(filedir+"/JERProviderPlots.root")
#mmcTool.SetAlgorithmVersion(1)
#mmcTool.SetNiterFit1(20)
#mmcTool.SetNiterFit2(30)
#mmcTool.SetNiterFit3(10)
#mmcTool.SetNsigmaMETscan(3.0)
#mmcTool.SetApplyMassScale(0)


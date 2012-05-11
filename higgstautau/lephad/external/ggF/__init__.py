#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Compile and load the egamma EnergyRescaler correction tool
print "    -----Loading the ggFReweighting tool (and compiling if necessary)..."
compileC(filedir + "/ggFReweighting.cxx")
from ROOT import AnalysisFramework
ggFReweighting = AnalysisFramework.External.ggFReweighting

ggF_tool = {}
ggF_tool[100] = ggFReweighting("PowHeg", 100, "Mean", filedir + '/', "mc11")
ggF_tool[105] = ggFReweighting("PowHeg", 105, "Mean", filedir + '/', "mc11")
ggF_tool[110] = ggFReweighting("PowHeg", 110, "Mean", filedir + '/', "mc11")
ggF_tool[115] = ggFReweighting("PowHeg", 115, "Mean", filedir + '/', "mc11")
ggF_tool[120] = ggFReweighting("PowHeg", 120, "Mean", filedir + '/', "mc11")
ggF_tool[125] = ggFReweighting("PowHeg", 125, "Mean", filedir + '/', "mc11")
ggF_tool[130] = ggFReweighting("PowHeg", 130, "Mean", filedir + '/', "mc11")
ggF_tool[135] = ggFReweighting("PowHeg", 135, "Mean", filedir + '/', "mc11")
ggF_tool[140] = ggFReweighting("PowHeg", 140, "Mean", filedir + '/', "mc11")
ggF_tool[145] = ggFReweighting("PowHeg", 145, "Mean", filedir + '/', "mc11")
ggF_tool[150] = ggFReweighting("PowHeg", 150, "Mean", filedir + '/', "mc11")

def getggFTool(higgsMass):
    if higgsMass in ggF_tool:
        return ggF_tool[higgsMass]
    else:
        print "@@@@@ Higgs mass cannot be determined for ggF tool...using 100 GeV for now"
        print "      This message is normal for non-Higgs samples."
        return ggF_tool[100]

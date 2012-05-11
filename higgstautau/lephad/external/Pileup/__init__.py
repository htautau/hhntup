#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Load the pile-up correction tool
print "    -----Loading the pile-up reweighting tool (and compiling if necessary)..."  
compileC(filedir+"/TPileupReweighting.cxx")

from ROOT import AnalysisFramework

pileup_mc11a = AnalysisFramework.External.Root.TPileupReweighting("AdvancedPURT_mc11a")
pileup_mc11a.SetDataScaleFactors(0.83)
pileup_mc11a.AddConfigFile(filedir + "/mc11a_defaults.root")
pileup_mc11a.AddLumiCalcFile(filedir + "/ilumicalc_histograms_None_178044-191933_slimmed.root")
pileup_mc11a.SetUnrepresentedDataAction(2)
pileup_mc11a.SetDefaultChannel(0)
#pileup_mc11a.SetDataScaleFactors(0.83)
pileup_mc11a.initialize()

pileup_mc11b = AnalysisFramework.External.Root.TPileupReweighting("AdvancedPURT_mc11b")
pileup_mc11b.SetDataScaleFactors(0.83)
pileup_mc11b.AddConfigFile(filedir + "/mc11b_defaults.root")
pileup_mc11b.AddLumiCalcFile(filedir + "/ilumicalc_histograms_None_178044-191933_slimmed.root")
pileup_mc11b.SetUnrepresentedDataAction(2)
pileup_mc11b.SetDefaultChannel(0)
#pileup_mc11b.SetDataScaleFactors(0.83)
pileup_mc11b.initialize()

pileup_mc11c = AnalysisFramework.External.Root.TPileupReweighting("AdvancedPURT_mc11b")
pileup_mc11c.SetDataScaleFactors(1.)
pileup_mc11c.AddConfigFile(filedir + "/mc11b_defaults.root")
pileup_mc11c.AddLumiCalcFile(filedir + "/ilumicalc_histograms_None_178044-191933_slimmed.root")
pileup_mc11c.SetUnrepresentedDataAction(2)
pileup_mc11c.SetDefaultChannel(0)
pileup_mc11c.initialize()


def getPileupTool(a, b, c):
    if a:
        return pileup_mc11a
    elif b:
        return pileup_mc11b
    elif c:
        return pileup_mc11c
    else:
        print "@@@@@ MC production cannot be determined for Pileup tool...using mc11c for now"
        return pileup_mc11c

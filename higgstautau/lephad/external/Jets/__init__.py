#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Load the JER tool
print "    -----Loading the JER tool..."
compileC(filedir + "/JERProvider.cxx")
from ROOT import AnalysisFramework
jer_tool = AnalysisFramework.External.JERProvider("AntiKt4LCTopoJES", "Truth", filedir+"/JERProviderPlots.root")
jer_tool.init()

# Load the JetUncertainties tool
print "    -----Loading the JetUncertainties tool..."
compileC(filedir + "/JESUncertaintyProvider.cxx")
compileC(filedir + "/MultijetJESUncertaintyProvider.cxx")
jes_tool = AnalysisFramework.External.MultijetJESUncertaintyProvider("AntiKt4LCTopoJets", filedir+"/UnknownFlavourComp.root", filedir+"/MJESUncertainty.root", filedir+"/JESUncertainty.root")
jes_tool.includeFlavorComposition(False)
jes_tool.init()

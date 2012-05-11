#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Load the GRL tool
print "    -----Loading the GoodRunsLists tool..."
compileC(filedir + "/TMsgLogger.cxx")
compileC(filedir + "/TLumiBlockRange.cxx")
compileC(filedir + "/TGoodRun.cxx")
compileC(filedir + "/TGoodRunsList.cxx")
compileC(filedir + "/TGRLCollection.cxx")
compileC(filedir + "/StrUtil.cxx")
compileC(filedir + "/TGoodRunsListReader.cxx")

from ROOT import AnalysisFramework
grlR = AnalysisFramework.External.TGoodRunsListReader()
grlR.SetXMLFile(filedir+'/data11_7TeV.periodAllYear_DetStatus-v36-pro10_CoolRunQuery-00-04-08_Higgs_tautau_lh.xml')
grlR.Interpret()
grl_tool = grlR.GetMergedGoodRunsList()
grl_tool.Summary(False)

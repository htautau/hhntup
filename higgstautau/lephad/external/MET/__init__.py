#Author: Kong Guan Tan <kong.guan.tan@cern.ch>

from ..Core import compileC
import os

filedir = os.path.dirname(os.path.abspath(__file__))

# Load the MET utility
print "    -----Loading the MET utility..."
compileC(filedir + "/METUtility.cxx")
from ROOT import AnalysisFramework
METUtility = AnalysisFramework.External.METUtility

def SetMETTerms(met_utility, met_object):
    met_utility.setMETTerm("RefJet",met_object.RefJet_etx(),met_object.RefJet_ety(),met_object.RefJet_sumet())
    met_utility.setMETTerm("RefEle",met_object.RefEle_etx(),met_object.RefEle_ety(),met_object.RefEle_sumet())
    met_utility.setMETTerm("RefMuon",met_object.RefMuon_etx(),met_object.RefMuon_ety(),met_object.RefMuon_sumet())
    met_utility.setMETTerm("RefTau",met_object.RefTau_etx(),met_object.RefTau_ety(),met_object.RefTau_sumet())
    met_utility.setMETTerm("CellOutEflow",met_object.CellOut_etx(),met_object.CellOut_ety(),met_object.CellOut_sumet())
    met_utility.setMETTerm("MuonTotal",met_object.Muon_Total_etx(),met_object.Muon_Total_ety(),met_object.Muon_Total_sumet())
    met_utility.setMETTerm("SoftJets",met_object.SoftJets_etx(),met_object.SoftJets_ety(),met_object.SoftJets_sumet())
    met_utility.setMETTerm("RefGamma",met_object.RefGamma_etx(),met_object.RefGamma_ety(),met_object.RefGamma_sumet())

def SetTauParameters(met_utility, met_object, tau_object, tau_num=0):
    met_utility.setTauParameters(tau_object.pt()[0],tau_object.eta()[0],tau_object.phi()[0],met_object.tau_MET_wet(),met_object.tau_MET_wpx(),met_object.tau_MET_wpy(),met_object.tau_MET_statusWord())

def DoMETSystematics(met_utility, met_object, option):
    SetMETTerms(met_utility, met_object)
    new_met = met_utility.getMissingET("RefFinal", option)
    met_object.RefFinal_etx( new_met.etx() )
    met_object.RefFinal_ety( new_met.ety() )
    met_object.RefFinal_phi( new_met.phi() )
    met_object.RefFinal_et( new_met.et() )
    met_object.RefFinal_sumet( new_met.sumet() )

met_utility = METUtility()

#met_utility.setObjectEnergyUncertainties("taus",tesUp,tesDown)
#met_utility.setObjectEnergyUncertainties("taus",tesUp,tesDown)

#met_utility.getMissingET("RefFinal","AllClustersUp").etx() # or Down
#met_utility.getMissingET("RefFinal","AllClustersUp").ety() # or Down

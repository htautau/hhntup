from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from math import *

"""
https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MSSMHiggsToTauTauToHH2011Summer
https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MSSMHiggsToTauTauToLH2011Summer
"""

class TauElectronVeto(EventFilter):
    
    def passes(self, event):

        event.taus.select(lambda tau: tau.BDTEleScore > 0.51)
        return len(event.taus) > 0


class TauMuonVeto(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.muonVeto == 0)
        return len(event.taus) > 0


class TwoGoodTaus(EventFilter):

    def passes(self, event):

        # Only consider taus with at least a calo seed and which have at least one track
        # in a kinematic region passing the tau muon veto and charge req.
        event.taus.select(lambda tau: tau.author != 2 and tau.seedCalo_numTrack > 0 and
                                      tau.pt > 20*GeV and
                                      abs(tau.charge) == 1)
        return len(event.taus) > 1


class PriVertex(EventFilter):

    def passes(self, event):

        return any(ifilter(lambda vertex: vertex.nTracks >= 4 and vertex.type == 1, event.vertices))


class OrigTrigger(EventFilter):

    def passes(self, event):
        
        return event.EF_tau29_medium1_tau20_medium1


class Trigger(EventFilter):

    def passes(self, event):
        
        return event.EF_tau29_medium1_tau20_medium1 or \
               event.EF_e20_medium or \
               event.EF_e60_loose or \
               event.EF_xe60_noMu


class TriggerNoXE(EventFilter):

    def passes(self, event):
        
        return event.EF_tau29_medium1_tau20_medium1 or \
               event.EF_e20_medium or \
               event.EF_e60_loose


class MET(EventFilter):

    def passes(self, event):
        
        METx = event.MET_LocHadTopo_etx + event.MET_MuonBoy_etx - event.MET_RefMuon_Track_etx
        METy = event.MET_LocHadTopo_ety + event.MET_MuonBoy_ety - event.MET_RefMuon_Track_ety
        MET = math.sqrt(METx**2 + METy**2)
        return MET > 25*GeV


class JetCleaningLoose(EventFilter):
    """
    https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HowToCleanJets#Bad_jets_rel16_data
    """
    def passes(self, event):

        for jet in event.jets:
            if jet.pt <= 20*GeV: continue
            # HEC spike
            if jet.hecf>0.5 and abs(jet.HECQuality)>0.5: return False
            if abs(jet.NegativeE) > 60*GeV: return False
            # EM coherent noise
            if jet.emfrac > 0.95 and abs(jet.LArQuality)>0.8 and abs(jet.emscale_eta) < 2.8: return False
            # Cosmics, beam background
            if abs(jet.Timing) > 25.: return False
            if jet.emfrac < 0.05 and jet.sumPtTrk/jet.pt < 0.05 and abs(jet.emscale_eta) < 2: return False
            if jet.emfrac < 0.05 and abs(jet.emscale_eta) >= 2: return False
            if jet.fracSamplingMax>0.99 and abs(jet.emscale_eta)<2.0: return False
         
        return True


class JetCleaningMedium(EventFilter):
    """
    https://twiki.cern.ch/twiki/bin/view/AtlasProtected/HowToCleanJets#Bad_jets_rel16_data
    """
    def passes(self, event):

        for jet in event.jets:
            if jet.pt <= 20*GeV: continue
            # HEC spike
            if jet.hecf > 1. - abs(jet.HECQuality): return False
            # EM coherent noise
            if jet.emfrac > 0.9 and abs(jet.LArQuality)>0.8 and abs(jet.emscale_eta) < 2.8: return False
            # Cosmics, beam background
            if abs(jet.Timing) > 10.: return False
            if jet.emfrac < 0.05 and jet.sumPtTrk/jet.pt < 0.1 and abs(jet.emscale_eta) < 2: return False
            if jet.emfrac > 0.95 and jet.sumPtTrk/jet.pt < 0.05 and abs(jet.emscale_eta) < 2: return False
         
        return True


class LArError(EventFilter):

    def passes(self, event):

        return event.larError == 0


class LArHole(EventFilter):

    def passes(self, event):

        return not any(ifilter(lambda jet: jet.pt > 40*GeV and (-0.1 < jet.emscale_eta < 1.5) and (-0.9 < jet.emscale_phi < -0.5), event.jets))


class JetCrackVeto(EventFilter):

    def passes(self, event):

        for jet in event.jets:
            if jet.pt <= 20*GeV: continue
            if 1.3 < abs(jet.emscale_eta) < 1.7: return False
        return True


class ElectronVeto(EventFilter):

    def passes(self, event):

       for el in event.electrons:
           pt = el.cl_E / cosh(el.tracketa)
           if pt <= 15*GeV: continue
           if not (abs(el.cl_eta) < 1.37 or 1.52 < abs(el.cl_eta) and abs(el.cl_eta) < 2.47): continue
           if not (el.OQ & 1446) == 0: continue
           if el.author not in (1, 3): continue
           if el.medium != 1: continue
           return False

       return True


class MuonVeto(EventFilter):

    def passes(self, event):

       for muon in event.muons:
           if muon.pt <= 10*GeV: continue
           if muon.author != 6: continue # STACO
           if muon.isCombinedMuon != 1: continue
           if abs(muon.eta) >= 2.5: continue
           if muon.tight != 1: continue
           return False
       
       return True

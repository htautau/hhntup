from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV

try:
    from atlasutils.jets import cleaning as jetcleaning
except ImportError:
    class jetcleaning(object):
        LOOSER = 0
    print "Could not import atlastools.jets.cleaning"

from atlastools import datasets
from math import *

"""
See main documentation here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter
"""

MIN_TAUS = 1


class TauElectronVeto(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.EleBDTLoose == 0)
        return len(event.taus) >= MIN_TAUS


class TauMuonVeto(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.muonVeto == 0)
        return len(event.taus) >= MIN_TAUS


class TauHasTrack(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.numTrack > 0)
        return len(event.taus) >= MIN_TAUS


class TauAuthor(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.author != 2)
        return len(event.taus) >= MIN_TAUS


class TauPT(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.pt > 20*GeV)
        return len(event.taus) >= MIN_TAUS


class TauEta(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: abs(tau.eta) < 2.1)
        return len(event.taus) >= MIN_TAUS


class TauJVF(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: True if abs(tau.track_eta[0]) > 2.1 else tau.jet_jvtxf > .5)
        return len(event.taus) >= MIN_TAUS


class Tau1Track3Track(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.numTrack in (1, 3))
        return len(event.taus) >= MIN_TAUS


class TauCharge(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: abs(tau.charge) == 1)
        return len(event.taus) >= MIN_TAUS


class TauLoose(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.JetBDTSigLoose)
        return len(event.taus) >= MIN_TAUS


class TauCrack(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: not (1.37 < abs(tau.track_eta[0]) < 1.52))
        return len(event.taus) >= MIN_TAUS


class TauLArHole(EventFilter):

    def passes(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        for tau in event.taus:
            if tau.track_n > 0:
                eta = tau.track_eta[0]
                phi = tau.track_phi[0]

                if (-0.1 < eta < 1.55) and (-0.9 < phi < -0.5):
                    return False

        return True

        # event.taus.select(lambda tau:
        #         not ((tau.track_n > 0)
        #              and (-0.1 < tau.track_eta[0] < 1.55)
        #              and (-0.9 < tau.track_phi[0] < -0.5)))
        # return len(event.taus) >= 1


class ElectronLArHole(EventFilter):

    def passes(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        for electron in event.electrons:
            eta = electron.cl_eta
            phi = electron.cl_phi

            if (-0.1 < eta < 1.55) and (-0.888 < phi < -0.492):
                return False

        return True


class muTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-H  (177986-184169) : EF_mu18_MG
        period I-M  (185353-191933) : EF_mu18_MG_medium
        """
        try:
            if 177986 <= event.RunNumber <= 184169: # Period B-H 1237.26 pb-1
                return event.EF_mu18_MG
            elif 185353 <= event.RunNumber <= 191933: # Period I-M 3982.9 pb-1
                return event.EF_mu18_MG_medium
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


data_triggers = [
    'EF_mu18_MG',
    'EF_mu18_MG_medium',
    'EF_e20_medium',
    'EF_e22_medium',
    'EF_e22vh_medium1'
]

mc_triggers = [
    'EF_mu18_MG',
    'EF_mu18_MG_medium',
    'EF_e20_medium',
    'EF_e22_medium',
    'EF_e22vh_medium1'
]


class muMCTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-H  (177986-187815) : EF_mu18_MG
        period I-M  (188902-191933) : EF_mu18_MG_medium
        """
        try:
            if 177986 <= event.RunNumber <= 187815: # Period B-H 1237.26 pb-1
                return event.EF_mu18_MG
            elif 188902 <= event.RunNumber <= 191933: # Period I-M 3982.9 pb-1
                return event.EF_mu18_MG_medium
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


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
           if not (abs(el.cl_eta) < 1.37 or 1.52 < abs(el.cl_eta) < 2.47): continue
           if not (el.OQ & 1446) == 0: continue
           if el.author not in (1, 3): continue
           if el.mediumPP != 1: continue
           if not abs(el.charge) == 1: continue
           return False

       return True


def muon_has_good_track(muon):

    blayer = (muon.expectBLayerHit == 0.0) or (muon.nBLHits > 0)
    pix = (muon.nPixHits + muon.nPixelDeadSensors) > 1
    sct = (muon.nSCTHits + muon.nSCTDeadSensors) > 5
    holes = (muon.nPixHoles + muon.nSCTHoles) < 3
    trt = False

    if abs(muon.eta) < 1.9:
        trt = (muon.nTRTHits + muon.nTRTOutliers) > 5 and \
              muon.nTRTOutliers < 0.9*(muon.nTRTHits + muon.nTRTOutliers)
    else:
        trt = (muon.nTRTHits + muon.nTRTOutliers) <= 5 or \
              (muon.nTRTOutliers) < 0.9 * ((muon.nTRTHits) + (muon.nTRTOutliers))

    return blayer and pix and sct and holes and trt


def tau_preselection(tau):
    """ Does the complete tau preselection """

    if not (tau.pt > 20*GeV) : return False
    if not (tau.numTrack == 1 or tau.numTrack == 3) : return False
    if not (tau.JetBDTSigMedium == 1) : return False
    if not (abs(tau.eta) < 2.5) : return False
    if not (tau.author != 2) : return False
    if not (abs(tau.charge) == 1) : return False
    if not (tau.EleBDTMedium == 0) : return False
    if not (tau.muonVeto == 0) : return False

    return True


def tau_selection(tau):
    """ Finalizes the tau selection """

    return True


def muon_preselection(mu):
    """ Does the complete muon preselection """

    if not (mu.pt > 10*GeV) : return False
    if not (muon_has_good_track(mu)) : return False
    if not (abs(mu.eta) < 2.5) : return False
    if not (mu.loose) : return False

    return True


def muon_selection(mu):
    """ Finalizes the muon selection """

    if not (mu.isCombinedMuon) : return False
    if not (mu.EFMG_matched) : return False

    # Find out if the muon is matched to a good trigger muon
    nMatches = mu.EFMG_n
    goodMatch = False
    for i in range(0, nMatches):
        if mu.EFMG_pt[i] > 18000: goodMatch = True

    return goodMatch




def electron_preselection(e):
    """ Does the complete electron preselection """

    cl_eta = e.cl_eta
    trk_eta = e.tracketa
    cl_Et = e.cl_E/cosh(trk_eta)
    
    if not (cl_Et > 15*GeV) : return False
    if not (abs(cl_eta) < 2.47) : return False
    if (1.37 < abs(cl_eta) < 1.52) : return False
    if not (e.author == 1 or e.author == 3) : return False
    if not (e.mediumPP) : return False

    return True


def electron_selection(e):
    """ Finalizes the electron selection"""

    if not (cl_Et > 20*GeV) : return False
    if not (e.tightPP) : return False

    return True


def jet_preselection(jet):
    """ Does the complete jet preselection """

    if not (jet.pt > 20*GeV) : return False

    return True


def jet_selection(jet):
    """ Finalizes the jet selection """

    if not (jet.pt > 25*GeV) : return False
    if not (abs(jet.eta) < 4.5) : return False
    if (abs(jet.eta) < 2.4):
        if not (jet.jvtxf > 0.75) : return False

    return True


def vertex_selection(vxp):
    """ Does the full primary and pileup vertex selection """

    return (vxp.type == 1 and vxp.nTracks >= 4) or (vxp.type == 3 and vxp.nTracks >= 2)


class MuonPreSelection(EventFilter):
    """Selects muons of good quality"""

    def passes(self, event):

        event.muons.select(lambda muon : muon_preselection(muon))
        return True


class MuonSelection(EventFilter):
    """Selects muons of good quality"""

    def passes(self, event):

        event.muons.select(lambda muon : muon_selection(muon))
        return len(event.muons) > 0


class ElectronPreSelection(EventFilter):
    """Selects electrons of good quality"""

    def passes(self, event):

        event.electrons.select(lambda electron : electron_preselection(electron))
        return True


class ElectronSelection(EventFilter):
    """Selects electrons of good quality"""

    def passes(self, event):

        event.electrons.select(lambda electron : electron_selection(electron))
        return len(event.electrons) > 0


class TauPreSelection(EventFilter):
    """Selects taus of good quality"""
    
    def passes(self, event):

        event.taus.select(lambda tau : tau_preselection(tau))
        return len(event.taus) > 0


class TauSelection(EventFilter):
    """Selects taus of good quality"""
    
    def passes(self, event):

        event.taus.select(lambda tau : tau_selection(tau))
        return len(event.taus) == 1


class JetPreSelection(EventFilter):
    """Selects jets of good quality, keep event in any case"""

    def passes(self, event):

        event.jets.select(lambda jet : jet_preselection(jet))
        return True


class JetSelection(EventFilter):
    """Selects jets of good quality, keep event in any case"""

    def passes(self, event):

        event.jets.select(lambda jet : jet_selection(jet))
        return True


class JetOverlapRemoval(EventFilter):
    """Muons > Electrons > Taus > Jets"""

    def passes(self, event):
        """ Remove jets matching muons and electrons """
        event.jets.select(lambda jet: not any([muon for muon in event.muons if (utils.dR(muon.eta, muon.phi, jet.eta, jet.phi) < 0.2)]))
        event.jets.select(lambda jet: not any([electron for electron in event.electrons if (utils.dR(electron.eta, electron.phi, jet.eta, jet.phi) < 0.2)]))
        return True


class LeptonOverlapRemoval(EventFilter):
    """ Muons > Electrons """

    def passes(self, event):
        """ Remove electrons matching muons """

        event.electrons.select(lambda electron: not any([muon for muon in event.muons if (utils.dR(muon.eta, muon.phi, electron.eta, electron.phi) < 0.2)]))

        return True


class FinalOverlapRemoval(EventFilter):
    """Muons > Electrons > Taus > Jets"""

    def passes(self, event):

        event.taus.select(lambda tau: not any([electron for electron in event.electrons if (utils.dR(electron.eta, electron.phi, tau.eta, tau.phi) < 0.2)]))
        event.taus.select(lambda tau: not any([muon for muon in event.muons if (utils.dR(muon.eta, muon.phi, tau.eta, tau.phi) < 0.2)]))
        event.jets.select(lambda jet: not any([tau for tau in event.taus if (utils.dR(tau.eta, tau.phi, jet.eta, jet.phi) < 0.2)]))

        return len(event.taus) == 1
            


class DileptonVeto(EventFilter):
    """Keep events with one muon only"""

    def passes(self, event):
        return len(event.muons) + len(event.electrons) == 1


class HasTau(EventFilter):
    """Keep events with at least one good tau"""

    def passes(self, event):
        return len(event.taus) > 0
        

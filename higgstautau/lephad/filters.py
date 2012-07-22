from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from rootpy.math.physics.vector import LorentzVector

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


class SetElectronsFourVector(EventFilter):

    def passes(self, event):
        
        for el in event.electrons:
            if ((el.nSCTHits + el.nPixHits) < 4):
                eta = el.cl_eta
                phi = el.cl_phi
                et  = el.cl_pt
            else:
                eta = el.tracketa
                phi = el.trackphi
                et  = el.cl_E / cosh(el.tracketa)
            vect = LorentzVector()
            vect.SetPtEtaPhiE(et, eta, phi, el.cl_E)
            setattr(el, 'fourvect', vect)

        return True


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



class JetCrackVeto(EventFilter):

    def passes(self, event):

        for jet in event.jets:
            if jet.pt <= 20*GeV: continue
            if 1.3 < abs(jet.emscale_eta) < 1.7: return False
        return True



############################################################
# TRIGGERS
############################################################

class noTriggers(EventFilter):
    """ Place holder for no triggers for embedding samples"""

    def passes(self, event):
        return True

#--------------------------------------------
# Muon Data Triggers
#--------------------------------------------

# Muon Single Lepton Triggers
class muSLTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-I  (177986-186493) : EF_mu18_MG
        period J-M  (186516-191933) : EF_mu18_MG_medium
        """
        try:
            if 177986 <= event.RunNumber <= 186493:
                return event.EF_mu18_MG
            elif 186516 <= event.RunNumber <= 191933:
                return event.EF_mu18_MG_medium
            elif 200000 <= event.RunNumber:
                return True
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Muon combined triggers
class muLTTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers#Tau_Lepton_Triggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-J  (177986-186755) : EF_tau16_loose_mu15
        period K-M  (186873-191933) : EF_tau20_medium_mu15
        """
        try:
            if 177986 <= event.RunNumber <= 186755:
                return event.EF_tau16_loose_mu15
            elif 186873 <= event.RunNumber <= 191933:
                return event.EF_tau20_medium_mu15
            elif 200000 <= event.RunNumber:
                return True
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Muon any trigger
class AnyMuTriggers(EventFilter):

    def passes(self, event):

        TriggerList = ['EF_mu18',
                       'EF_mu18_MG',
                       'EF_mu18_medium',
                       'EF_mu18_MG_medium',
                       'EF_tau16_loose_mu15',
                       'EF_tau20_medium_mu15']

        TriggersToOR = 0

        for trig in TriggerList:
            try:
                TriggersToOR += getattr(event, trig)
            except AttributeError:
                pass
                
        if TriggersToOR > 0: return True
        else: return False


# Muon SLT or LLT
class AllMuTriggers(EventFilter):

    def passes(self, event):
        SLT = muSLTriggers()
        LLT = muLTTriggers()

        return SLT.passes(event) or LLT.passes(event)



#--------------------------------------------
# Electron Data Triggers
#--------------------------------------------

# Electron single lepton triggers
class eSLTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-J  (177986-186755) : EF_e20_medium
        period K    (186873-187815) : EF_e22_medium
        period L-M  (188902-191933) : EF_e22vh_medium1
        """
        try:
            if 177986 <= event.RunNumber <= 186755:
                return event.EF_e20_medium
            elif 186873 <= event.RunNumber <= 187815:
                return event.EF_e22_medium
            elif 188902 <= event.RunNumber <= 191933:
                return event.EF_e22vh_medium1
            elif 200000 <= event.RunNumber:
                return True
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Electron combined triggers
class eLTTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers#Tau_Lepton_Triggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-J  (177986-186755) : EF_tau16_e15_medium
        period K    (186873-187815) : EF_tau20_medium_e15_medium
        period L-M  (188902-191933) : EF_tau20_medium_e15vh_medium
        """
        try:
            if 177986 <= event.RunNumber <= 186755:
                return event.EF_tau16_e15_medium
            elif 186873 <= event.RunNumber <= 187815:
                return event.EF_tau20_medium_e15_medium
            elif 188902 <= event.RunNumber <= 191933:
                return event.EF_tau20_medium_e15vh_medium
            elif 200000 <= event.RunNumber:
                return True
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Electron any trigger
class AnyETriggers(EventFilter):

    def passes(self, event):

        TriggerList = ['EF_e20_medium',
                       'EF_e22_medium',
                       'EF_e22vh_medium1',
                       'EF_tau16_loose_e15_medium',
                       'EF_tau20_medium_e15_medium',
                       'EF_tau20_medium_e15vh_medium']

        TriggersToOR = 0

        for trig in TriggerList:
            try:
                TriggersToOR += getattr(event, trig)
            except AttributeError:
                pass
                
        if TriggersToOR: return True
        else: return False


# Electron SLT or LLT
class AllETriggers(EventFilter):

    def passes(self, event):
        SLT = eSLTriggers()
        LLT = eLTTriggers()

        return SLT.passes(event) or LLT.passes(event)


#--------------------------------------------
# Muon MC Triggers
#--------------------------------------------

            
# Muon single lepton triggers
class muMCSLTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-K  (177986-187815) : EF_mu18_MG
        period L-M  (188902-191933) : EF_mu18_MG_medium
        """
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_mu18_MG
            elif 188902 <= event.RunNumber:
                return event.EF_mu18_MG_medium
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Muon combined triggers
class muMCLTTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers#Tau_Lepton_Triggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-K  (177986-187815) : EF_tau16_loose_mu15
        period L-M  (188902-191933) : EF_tau20_medium_mu15
        """
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_tau16_loose_mu15
            elif 188902 <= event.RunNumber:
                return event.EF_tau20_medium_mu15
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


#--------------------------------------------
# Electron MC Triggers
#--------------------------------------------

# Electron single lepton triggers
class eMCSLTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-K  (177986-187815) : EF_e20_medium
        period L-M  (188902-191933) : EF_e22vh_medium1
        """
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_e20_medium
            elif 188902 <= event.RunNumber:
                return event.EF_e22vh_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Electron combined triggers
class eMCLTTriggers(EventFilter):

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers#Tau_Lepton_Triggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        period B-K  (177986-187815) : EF_tau16_loose_e15_medium
        period L-M  (188902-191933) : EF_tau20_medium_e15vh_medium
        """
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_tau16_loose_e15_medium
            elif 188902 <= event.RunNumber:
                return event.EF_tau20_medium_e15vh_medium
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


#--------------------------------------------
# Lepton MC Triggers
#--------------------------------------------

# Apply MC triggers
class AllMCTriggers(EventFilter):

    def passes(self, event):
        muonSLT = muMCSLTriggers()
        elecSLT = eMCSLTriggers()
        muonLTT = muMCLTTriggers()
        elecLTT = eMCLTTriggers()
        
        return muonSLT.passes(event) or elecSLT.passes(event) or muonLTT.passes(event) or elecLTT.passes(event)


# Use any trigger
class AnyMCTriggers(EventFilter):

    def passes(self, event):
        muonTrig = AnyMuTriggers()
        elecTrig = AnyETriggers()
        
        return muonTrig.passes(event) or elecTrig.passes(event)



############################################################
# TAU SELECTION
############################################################

def tau_skimselection(tau):
    """ Does the complete tau preselection """

    if not (tau.pt > 15*GeV) : return False
    if not (tau.numTrack == 1 or tau.numTrack == 3) : return False
    if not (tau.JetBDTSigLoose == 1) : return False
    if not (abs(tau.eta) < 2.5) : return False

    return True


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


############################################################
# MUON SELECTION
############################################################

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


def muon_skimselection(mu):
    """ Does the complete muon preselection """

    if not (mu.pt > 13*GeV) : return False
    if not (abs(mu.eta) < 3.0) : return False
    if not (mu.loose) : return False

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
    if not (mu.pt > 20*GeV) : return False

    return True


############################################################
# ELECTRON SELECTION
############################################################

def electron_skimselection(el):
    """ Does the complete electron preselection """

    el_cl_Et = getattr(el,'fourvect').Pt()

    if not (el_cl_Et > 15*GeV) : return False
    if not (abs(el.cl_eta) < 3.0) : return False
    if not (el.author == 1 or el.author == 3) : return False
    if not (el.mediumPP or el.tightPP) : return False

    return True


def electron_preselection(el):
    """ Does the complete electron preselection """

    el_cl_Et = getattr(el,'fourvect').Pt()

    if not (el_cl_Et > 15*GeV) : return False
    if not (abs(el.cl_eta) < 2.47) : return False
    if (1.37 < abs(el.cl_eta) < 1.52) : return False
    if not (el.author == 1 or el.author == 3) : return False
    if not (el.mediumPP) : return False

    return True


def electron_selection(el):
    """ Finalizes the electron selection"""

    el_cl_Et = getattr(el,'fourvect').Pt()

    if not (el_cl_Et > 25*GeV) : return False
    if not (el.tightPP) : return False

    return True


############################################################
# JET SELECTION
############################################################

def jet_preselection(jet):
    """ Does the complete jet preselection """

    if not (jet.pt > 20*GeV) : return False

    return True


def jet_selection(jet):
    """ Finalizes the jet selection """

    if not (jet.pt > 25*GeV) : return False
        
    #Protection against bunny ear jets
    if (2.5 < abs(jet.eta) < 3.5):
        if not (jet.pt > 30*GeV) : return False

    if not (abs(jet.eta) < 4.5) : return False
    if (abs(jet.eta) < 2.4):
        if not (jet.jvtxf > 0.75) : return False

    return True


############################################################
# VERTEX SELECTION
############################################################

def vertex_selection(vxp):
    """ Does the full primary and pileup vertex selection """

    return (vxp.type == 1 and vxp.nTracks >= 4) or (vxp.type == 3 and vxp.nTracks >= 2)



############################################################
# OBJECT ANALYSIS FILTERS
############################################################

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


############################################################
# OBJECT OVERLAP TOOLS
############################################################

def OverlapCheck(event, DoMuonCheck = False, DoElectronCheck = False):
    """
    Check for overlap between object, so that not a single object fires both
    the tau and the lepton requirements
    """
    if DoElectronCheck:
        for el in event.electrons:
            for tau in event.taus:
                if utils.dR(getattr(el,'fourvect').Eta(), getattr(el,'fourvect').Phi(), tau.eta, tau.phi) > 0.2:
                    return True
                
    if DoMuonCheck:
        for mu in event.muons:
            for tau in event.taus:
                if utils.dR(mu.eta, mu.phi, tau.eta, tau.phi) > 0.2:
                    return True

    return False

    
class JetOverlapRemoval(EventFilter):
    """Muons > Electrons > Taus > Jets"""

    def passes(self, event):
        """ Remove jets matching muons and electrons """
        event.jets.select(lambda jet: not any([mu for mu in event.muons if (utils.dR(mu.eta, mu.phi, jet.eta, jet.phi) < 0.2)]))
        event.jets.select(lambda jet: not any([el for el in event.electrons if (utils.dR(getattr(el,'fourvect').Eta(), getattr(el,'fourvect').Phi(), jet.eta, jet.phi) < 0.2)]))
        return True


class LeptonOverlapRemoval(EventFilter):
    """ Muons > Electrons """

    def passes(self, event):
        """ Remove electrons matching muons """

        event.electrons.select(lambda el: not any([mu for mu in event.muons if (utils.dR(mu.eta, mu.phi, getattr(el,'fourvect').Eta(), getattr(el,'fourvect').Phi()) < 0.2)]))

        return True


class FinalOverlapRemoval(EventFilter):
    """Muons > Electrons > Taus > Jets"""

    def passes(self, event):

        event.taus.select(lambda tau: not any([el for el in event.electrons if (utils.dR(getattr(el,'fourvect').Eta(), getattr(el,'fourvect').Phi(), tau.eta, tau.phi) < 0.2)]))
        event.taus.select(lambda tau: not any([mu for mu in event.muons if (utils.dR(mu.eta, mu.phi, tau.eta, tau.phi) < 0.2)]))
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
        

#################################################
# Electron ID recalculation (p1130 tau D3PDs)
#################################################

from externaltools import egammaAnalysisUtils
from ROOT import isTightPlusPlus, isMediumPlusPlus, egammaMenu

class ElectronIDpatch(EventFilter):
    """
    Recalculates the electron ID based on electron variables
    """

    def passes(self, event):

        #iterate over all electron candidates
        for e in event.electrons:

            #Retrieve relevant branches
            eta          = float(e.etas2)
            eT           = e.cl_E/cosh(eta)
            rHad         = e.Ethad/eT
            rHad1        = e.Ethad1/eT
            Reta         = e.reta
            w2           = e.weta2
            f1           = e.f1
            f3           = e.f3
            wstot        = e.wstot
            DEmaxs1      = -1
            if (e.emaxs1 + e.Emax2) != 0:
                DEmaxs1      = (e.emaxs1 - e.Emax2)/(e.emaxs1 + e.Emax2)
            deltaEta     = e.deltaeta1
            deltaPhi     = e.deltaphi2
            eOverp       = e.cl_E*abs(e.trackqoverp)
            d0           = e.trackd0_physics
            TRratio      = e.TRTHighTOutliersRatio
            nTRT         = e.nTRTHits
            nTRTOutliers = e.nTRTOutliers
            nSi          = e.nSiHits
            nSiOutliers  = e.nSCTOutliers + e.nPixelOutliers
            nPix         = e.nPixHits
            nPixOutliers = e.nPixelOutliers
            nBlayer      = e.nBLHits
            nBlayerOutliers = e.nBLayerOutliers
            expectBlayer = bool(e.expectHitInBLayer)
            ConvBit      = e.isEM

            #Correct the mediumPP flag
            e.mediumPP  = isMediumPlusPlus(eta, eT, f3, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1, deltaEta, d0, TRratio, nTRT, nTRTOutliers,
                                           nSi, nSiOutliers, nPix, nPixOutliers, nBlayer, nBlayerOutliers, expectBlayer)


            #Correct the tightPP flag
            e.tightPP = isTightPlusPlus(eta, eT, f3, rHad, rHad1, Reta, w2, f1, wstot, DEmaxs1, deltaEta, d0, TRratio, nTRT, nTRTOutliers,
                                        nSi, nSiOutliers, nPix, nPixOutliers, nBlayer, nBlayerOutliers, expectBlayer, eOverp, deltaPhi, ConvBit)
            
        return True
    


#################################################
# Tau ID recalculation (p1130 tau D3PDs)
#################################################

class TauIDpatch(EventFilter):
    """
    Recalculates the tau ID
    """

    def __init__(self, graph, **kwargs):

        super(TauIDpatch, self).__init__(**kwargs)

        #Load TGraphs

        f = graph

        self.loose_1p   = f.Get('loose_1p')
        self.medium_1p  = f.Get('medium_1p')
        self.tight_1p   = f.Get('tight_1p')
        self.loose_3p   = f.Get('loose_3p')
        self.medium_3p  = f.Get('medium_3p')
        self.tight_3p   = f.Get('tight_3p')


    def passes(self, event):

        for tau in event.taus:
            pt = tau.pt

            cut_loose  = 0
            cut_medium = 0
            cut_tight  = 0
            
            if tau.numTrack <= 1:
                cut_loose  = self.loose_1p.Eval(pt)
                cut_medium = self.medium_1p.Eval(pt)
                cut_tight  = self.tight_1p.Eval(pt)
            else:
                cut_loose  = self.loose_3p.Eval(pt)
                cut_medium = self.medium_3p.Eval(pt)
                cut_tight  = self.tight_3p.Eval(pt)

            tau.JetBDTSigLoose  = (tau.BDTJetScore > cut_loose)
            tau.JetBDTSigMedium = (tau.BDTJetScore > cut_medium)
            tau.JetBDTSigTight  = (tau.BDTJetScore > cut_tight)
        
        return True

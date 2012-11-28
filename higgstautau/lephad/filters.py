from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from rootpy.math.physics.vector import LorentzVector
from atlastools import datasets
from math import *

"""
See main documentation here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter
"""

MIN_TAUS = 1


class PrepareInputTree(EventFilter):

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
            vect.SetPtEtaPhiM(et, eta, phi, 0.510999)
            setattr(el, 'fourvect', vect)

        ## Append the isLTT (is Lepton-Tau trigger event) flag
        event.isLTT = False

        ## Default value for the leptonType
        event.leptonType = 'mu'

        return True


class AcceptanceChallenge(EventFilter):

    def passes(self, event):

        EventNumbers = [30099, 30100, 30122, 30138, 30166, 30217, 30303, 30311, 30343, 30347]
        EventNumbers += [30095, 30125, 30134, 30154, 30167, 30182, 30219, 30227, 30248, 30266]

        if event.EventNumber not in EventNumbers:
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




############################################################
# TRIGGERS
############################################################

class noTriggers(EventFilter):
    """ Place holder for no triggers for embedding samples"""

    def __init__(self, year, **kwargs):

        self.year = year
        super(noTriggers, self).__init__(**kwargs)

    def passes(self, event):
        return True

#--------------------------------------------
# Muon Data Triggers
#--------------------------------------------

# Muon Single Lepton Triggers
class muSLTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(muSLTTriggers, self).__init__(**kwargs)

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-I  (177986-186493) : EF_mu18_MG
        period J-M  (186516-191933) : EF_mu18_MG_medium

        2012:
        period A-E (200804-210308) : EF_mu24i_tight
        """
        try:
            if 177986 <= event.RunNumber <= 186493:
                return event.EF_mu18_MG
            elif 186516 <= event.RunNumber <= 191933:
                return event.EF_mu18_MG_medium
            elif 200804 <= event.RunNumber <= 210308 and self.year == 2012:
                return event.EF_mu24i_tight
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Muon combined triggers
class muLTTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(muLTTriggers, self).__init__(**kwargs)
    
    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-J  (177986-186755) : EF_tau16_loose_mu15
        period K-M  (186873-191933) : EF_tau20_medium_mu15

        2012:
        period A-E (200804-210308) : EF_tau20_medium1_mu15
        """
        try:
            if 177986 <= event.RunNumber <= 186755:
                return event.EF_tau16_loose_mu15
            elif 186873 <= event.RunNumber <= 191933:
                return event.EF_tau20_medium_mu15
            elif 200804 <= event.RunNumber <= 210308 and self.year == 2012:
                return event.EF_tau20_medium1_mu15
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Muon any trigger
class AnyMuTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(AnyMuTrigger, self).__init__(**kwargs)

    def passes(self, event):

        TriggerList = []
        
        if self.year == 2011:
            TriggerList = [
                'EF_mu18',
                'EF_mu18_MG',
                'EF_mu18_medium',
                'EF_mu18_MG_medium',
                'EF_tau16_loose_mu15',
                'EF_tau20_medium_mu15'
                ]
            
        elif self.year == 2012:
            TriggerList = [
                'EF_mu24i_tight',
                'EF_tau20_medium1_mu15'
            ]

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

    def __init__(self, year, **kwargs):

        self.year = year
        super(AllMuTriggers, self).__init__(**kwargs)

    def passes(self, event):
        SLT = muSLTriggers(self.year)
        LLT = muLTTriggers(self.year)

        return SLT.passes(event) or LLT.passes(event)



#--------------------------------------------
# Electron Data Triggers
#--------------------------------------------

# Electron single lepton triggers
class eSLTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(eSLTTriggers, self).__init__(**kwargs)
    
    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-J  (177986-186755) : EF_e20_medium
        period K    (186873-187815) : EF_e22_medium
        period L-M  (188902-191933) : EF_e22vh_medium1

        2012:
        period A-E (200804-210308) : EF_e24vhi_tight
        """
        
        try:
            if 177986 <= event.RunNumber <= 186755:
                return event.EF_e20_medium
            elif 186873 <= event.RunNumber <= 187815:
                return event.EF_e22_medium
            elif 188902 <= event.RunNumber <= 191933:
                return event.EF_e22vh_medium1
            elif 200804 <= event.RunNumber <= 210308 and self.year == 2012:
                return event.EF_e24vhi_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Electron combined triggers
class eLTTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(eLTTriggers, self).__init__(**kwargs)
    
    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-J  (177986-186755) : EF_tau16_e15_medium
        period K    (186873-187815) : EF_tau20_medium_e15_medium
        period L-M  (188902-191933) : EF_tau20_medium_e15vh_medium

        2012:
        period A-E (200804-210308) : EF_tau20Ti_medium1_e18vh_medium1
        """
        
        try:
            if 177986 <= event.RunNumber <= 186755:
                return event.EF_tau16_loose_e15_medium
            elif 186873 <= event.RunNumber <= 187815:
                return event.EF_tau20_medium_e15_medium
            elif 188902 <= event.RunNumber <= 191933:
                return event.EF_tau20_medium_e15vh_medium
            elif 200804 <= event.RunNumber <= 210308 and year == 2012:
                return event.EF_tau20Ti_medium1_e18vh_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Electron any trigger
class AnyETriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(AnyETrigger, self).__init__(**kwargs)
    
    def passes(self, event):

        TriggerList = []
        
        if self.year == 2011:
            TriggerList = [
                'EF_e20_medium',
                'EF_e22_medium',
                'EF_e22vh_medium1',
                'EF_tau16_loose_e15_medium',
                'EF_tau20_medium_e15_medium',
                'EF_tau20_medium_e15vh_medium'
                ]

        elif self.year == 2012:
            TriggerList = [
                'EF_e24vhi_medium1',
                'EF_tau20Ti_medium1_e18vh_medium1'
                ]

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

    def __init__(self, year, **kwargs):

        self.year = year
        super(AllETriggers, self).__init__(**kwargs)

    def passes(self, event):
        SLT = eSLTriggers(self.year)
        LLT = eLTTriggers(self.year)

        return SLT.passes(event) or LLT.passes(event)


#--------------------------------------------
# Combined Data Triggers
#--------------------------------------------
class AllDataLTTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(AllDataLTTriggers, self).__init__(**kwargs)
    
    def passes(self, event):
        muSLT = muSLTriggers(self.year)
        muLTT = muLTTriggers(self.year)
        eSLT  = eSLTriggers(self.year)
        eLTT  = eLTTriggers(self.year)

        isSLT = muSLT.passes(event) or eSLT.passes(event)
        isLTT = muLTT.passes(event) or eLTT.passes(event)

        return isLTT and not isSLT


#--------------------------------------------
# Muon MC Triggers
#--------------------------------------------


# Muon single lepton triggers
class muMCSLTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(muMCSLTriggers, self).__init__(**kwargs)

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLH2012Winter#Data
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-K  (177986-187815) : EF_mu18_MG
        period L-M  (188902-191933) : EF_mu18_MG_medium

        2012:
        period A-E  (200804-210308) : EF_mu24i_tight
        """
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_mu18_MG
            elif 188902 <= event.RunNumber <= 195000:
                return event.EF_mu18_MG_medium
            elif self.year == 2012:
                return event.EF_mu24i_tight
            
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Muon combined triggers
class muMCLTTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(muMCLTTriggers, self).__init__(**kwargs)

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-K  (177986-187815) : EF_tau16_loose_mu15
        period L-M  (188902-191933) : EF_tau20_medium_mu15

        2012:
        period A-E  (200804-210308) : EF_tau20_medium1_mu15
        """
        
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_tau16_loose_mu15
            elif 188902 <= event.RunNumber <= 195000:
                return event.EF_tau20_medium_mu15
            elif self.year == 2012:
                return event.EF_tau20_medium1_mu15
            
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


#--------------------------------------------
# Electron MC Triggers
#--------------------------------------------

# Electron single lepton triggers
class eMCSLTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(eMCSLTriggers, self).__init__(**kwargs)

    def passes(self, event):
        """
         https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-K  (177986-187815) : EF_e20_medium
        period L-M  (188902-191933) : EF_e22vh_medium1

        2012:
        period A-E  (200804-210308) : EF_e24vhi_medium1
        """
        
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_e20_medium
            elif 188902 <= event.RunNumber <= 195000:
                return event.EF_e22vh_medium1
            elif self.year == 2012:
                return event.EF_e24vhi_medium1
            
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


# Electron combined triggers
class eMCLTTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(eMCLTTriggers, self).__init__(**kwargs)

    def passes(self, event):
        """
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToLHTriggers
        https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/DataPeriods

        2011:
        period B-K  (177986-187815) : EF_tau16_loose_e15_medium
        period L-M  (188902-191933) : EF_tau20_medium_e15vh_medium

        2012:
        period A-E  (200804-210308) : EF_tau20_medium1_e18vh_medium1
        """
        
        try:
            if 177986 <= event.RunNumber <= 187815:
                return event.EF_tau16_loose_e15_medium
            elif 188902 <= event.RunNumber <= 195000:
                return event.EF_tau20_medium_e15vh_medium
            elif self.year == 2012:
                return event.EF_tau20Ti_medium1_e18vh_medium1
            
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


#--------------------------------------------
# Lepton MC Triggers
#--------------------------------------------

# Apply MC triggers
class AllMCTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(AllMCTriggers, self).__init__(**kwargs)

    def passes(self, event):
        muonSLT = muMCSLTriggers(self.year)
        elecSLT = eMCSLTriggers(self.year)
        muonLTT = muMCLTTriggers(self.year)
        elecLTT = eMCLTTriggers(self.year)

        isSLT = muonSLT.passes(event) or elecSLT.passes(event)
        isLTT = muonLTT.passes(event) or elecLTT.passes(event)

        return isSLT or isLTT


# Use any trigger
class AnyMCTriggers(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(AnyMCTriggers, self).__init__(**kwargs)

    def passes(self, event):
        muonTrig = AnyMuTriggers(self.year)
        elecTrig = AnyETriggers(self.year)

        return muonTrig.passes(event) or elecTrig.passes(event)



############################################################
# TAU SELECTION
############################################################

def tau_skimselection(tau):
    """ Does the complete tau preselection """

    if not (tau.pt > 15*GeV) : return False
    if not (tau.numTrack == 1 or tau.numTrack == 3) : return False
    #if not (tau.JetBDTSigLoose == 1) : return False
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
    if tau.numTrack == 1:
        if not (tau.EleBDTMedium == 0) : return False
    if not (tau.muonVeto == 0) : return False

    return True


def tau_selection(tau):
    """ Finalizes the tau selection """

    return True


############################################################
# MUON SELECTION
############################################################

# use common method:
from ..filters import muon_has_good_track


def muon_skimselection(mu):
    """ Does the complete muon preselection """

    if not (mu.pt > 13*GeV) : return False
    if not (abs(mu.eta) < 3.0) : return False
    if not (mu.loose) : return False

    return True


def muon_overlap_selection(mu, year):
    """ Do a pre-preselection for overlap removal with taus and jets """

    if not (mu.pt > 4*GeV) : return False
    if not (muon_has_good_track(mu, year)) : return False
    if not (abs(mu.eta) < 2.5) : return False
    if not (mu.loose) : return False

    return True


def muon_preselection(mu, year):
    """ Does the complete muon preselection """

    if not (mu.pt > 10*GeV) : return False
    if not (muon_has_good_track(mu, year)) : return False
    if not (abs(mu.eta) < 2.5) : return False
    if not (mu.loose) : return False

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


def electron_preselection(el, year):
    """ Does the complete electron preselection """

    el_cl_Et = getattr(el,'fourvect').Pt()

    if not (el_cl_Et > 15*GeV) : return False
    if year == 2012:
        if not (el.OQ&1446) == 0: return False
    if not (abs(el.cl_eta) < 2.47) : return False
    if not (el.author == 1 or el.author == 3) : return False
    if not (el.loosePP) : return False

    return True

def electron_etaselection(el):
    """ Does the complete electron preselection """

    if (1.37 < abs(el.cl_eta) < 1.52) : return False

    return True



############################################################
# JET SELECTION
############################################################

def jet_preselection(jet):
    """ Does the complete jet preselection """

    if not (jet.pt > 20*GeV) : return False

    return True

############################################################
# OBJECT ANALYSIS FILTERS
############################################################

"""
NOTE: no need for lambda functions below, just pass the selection function
directly. Right now you pass a function that calls a function...
"""

class MuonOverlapSelection(EventFilter):
    """ Selects low pt muons of good quality for overlap removal with taus """

    def __init__(self, year, **kwargs):

        self.year = year
        super(MuonOverlapSelection, self).__init__(**kwargs)

    def passes(self, event):

        event.muons.select(lambda muon : muon_overlap_selection(muon, self.year))
        return True


class MuonPreSelection(EventFilter):
    """Selects muons of good quality"""

    def __init__(self, year, **kwargs):

        self.year = year
        super(MuonPreSelection, self).__init__(**kwargs)

    def passes(self, event):

        event.muons.select(lambda muon : muon_preselection(muon, self.year))
        return True


class ElectronPreSelection(EventFilter):
    """Selects electrons of good quality"""

    def __init__(self, year, **kwargs):

        self.year = year
        super(ElectronPreSelection, self).__init__(**kwargs)

    def passes(self, event):

        event.electrons.select(lambda electron : electron_preselection(electron, self.year))
        return True


class ElectronEtaSelection(EventFilter):
    """Selects electrons outside of the crack region"""

    def passes(self, event):

        event.electrons.select(lambda electron : electron_etaselection(electron))
        return True


class LeptonSelection(EventFilter):
    """ Selects the lepton, with all possible trigger/stream provenance """

    def __init__(self, datatype, stream, year, **kwargs):

        self.datatype = datatype
        self.stream = stream
        self.year = year
        super(LeptonSelection, self).__init__(**kwargs)

    def passes(self, event):

        # Get the lepton Pt
        Pt = None
        ID = None
        if event.leptonType == 'e':
            Pt = getattr(event.electrons[0],'fourvect').Pt()
            ID = event.electrons[0].tightPP
        if event.leptonType == 'mu':
            Pt = event.muons[0].pt
            ID = event.muons[0].isCombinedMuon

        mu_SLT_threshold = 0
        mu_LTT_threshold = 0
        el_SLT_threshold = 0
        el_LTT_threshold = 0

        if self.year == 2011:
            mu_SLT_threshold = 22*GeV
            mu_LTT_threshold = 17*GeV
            el_SLT_threshold = 25*GeV
            el_LTT_threshold = 17*GeV

        if self.year == 2012:
            mu_SLT_threshold = 26*GeV
            mu_LTT_threshold = 17*GeV
            el_SLT_threshold = 26*GeV
            el_LTT_threshold = 20*GeV

        if not ID: return False

        if self.datatype == datasets.EMBED:
                 
            if event.leptonType == 'e':
                if Pt > el_LTT_threshold:
                    if Pt <= el_SLT_threshold:
                        event.isLTT = True
                    return True
                else:
                    return False

            if event.leptonType == 'mu':
                if Pt > mu_LTT_threshold:
                    if Pt <= mu_SLT_threshold:
                        event.isLTT = True
                    return True
                else:
                    return False

            
        if self.datatype == datasets.MC:
                 
            if event.leptonType == 'e':
                if Pt > el_SLT_threshold:
                    trigger = eMCSLTriggers(self.year)
                    if trigger.passes(event):
                        return True
                if el_LTT_threshold < Pt <= el_SLT_threshold:
                    trigger = eMCLTTriggers(self.year)
                    if trigger.passes(event):
                        event.isLTT = True
                        return True
                return False

            if event.leptonType == 'mu':
                if Pt > mu_SLT_threshold:
                    trigger = muMCSLTriggers(self.year)
                    if trigger.passes(event):
                        return True
                if mu_LTT_threshold < Pt <= mu_SLT_threshold:
                    trigger = muMCLTTriggers(self.year)
                    if trigger.passes(event):
                        event.isLTT = True
                        return True
                return False

            
        if self.datatype == datasets.DATA:
            
            if self.stream == 'JetTauEtmiss':
                
                if event.leptonType == 'e':
                    SLT = eSLTriggers(self.year)
                    LTT = eLTTriggers(self.year)
                    if LTT.passes(event) and not SLT.passes(event):
                        if el_LTT_threshold < Pt <= el_SLT_threshold:
                            event.isLTT = True
                            return True
                    return False

                if event.leptonType == 'mu':
                    SLT = muSLTriggers(self.year)
                    LTT = muLTTriggers(self.year)
                    if LTT.passes(event) and not SLT.passes(event):
                        if mu_LTT_threshold < Pt <= mu_SLT_threshold:
                            event.isLTT = True
                            return True
                    return False

            if self.stream == 'Egamma' and event.leptonType == 'e':
                SLT = eSLTriggers(self.year)
                LTT = eLTTriggers(self.year)
                isSLT = SLT.passes(event)
                isLTT = LTT.passes(event)
                if isSLT and Pt > el_SLT_threshold:
                    return True
                if isLTT and isSLT and el_LTT_threshold < Pt <= el_SLT_threshold:
                    event.isLTT = True
                    return True
                return False

            if self.stream == 'Muons' and event.leptonType == 'mu':
                SLT = muSLTriggers(self.year)
                LTT = muLTTriggers(self.year)
                isSLT = SLT.passes(event)
                isLTT = LTT.passes(event)
                if isSLT and Pt > mu_SLT_threshold:
                    return True
                if isLTT and isSLT and mu_LTT_threshold < Pt <= mu_SLT_threshold:
                    event.isLTT = True
                    return True
                return False
                    

class TauPreSelection(EventFilter):
    """Selects taus of good quality"""

    def passes(self, event):

        event.taus.select(lambda tau : tau_preselection(tau))
        return len(event.taus) > 0


class TauSelection(EventFilter):
    """Selects taus of good quality"""

    def passes(self, event):

        event.taus.select(lambda tau : tau_selection(tau))
        if len(event.taus) == 1:
            if event.isLTT:
                return event.taus[0].pt > 25*GeV
            else:
                return True
        return False


class JetPreSelection(EventFilter):
    """Selects jets of good quality, keep event in any case"""

    def passes(self, event):

        event.jets.select(lambda jet : jet_preselection(jet))
        return True


# JetSelection is in the common filters
# putting import here so I don't break lephad code for now...
from ..filters import JetSelection


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
    """Muons, Electrons > Jets"""

    def passes(self, event):
        """ Remove jets matching muons and electrons """
        event.jets.select(lambda jet: not any([mu for mu in event.muons if (utils.dR(mu.eta, mu.phi, jet.eta, jet.phi) < 0.2)]))
        event.jets.select(lambda jet: not any([el for el in event.electrons if (utils.dR(getattr(el,'fourvect').Eta(), getattr(el,'fourvect').Phi(), jet.eta, jet.phi) < 0.2)]))
        return True


class TauMuonOverlapRemoval(EventFilter):
    """Muons, Electrons > Jets"""

    def passes(self, event):
        """ Remove jets matching muons and electrons """
        event.taus.select(lambda tau: not any([mu for mu in event.muons if (utils.dR(mu.eta, mu.phi, tau.eta, tau.phi) < 0.2)]))
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
        event.jets.select(lambda jet: not any([tau for tau in event.taus if (utils.dR(tau.eta, tau.phi, jet.eta, jet.phi) < 0.2)]))

        return len(event.taus) == 1



## Other Event selection

class DileptonVeto(EventFilter):
    """Keep events with one muon only"""

    def passes(self, event):

        nMu = len(event.muons)
        nE  = len(event.electrons)

        if nE + nMu == 1:
            if nE > 0:
                event.leptonType = 'e'
            if nMu > 0:
                event.leptonType = 'mu'

            return True
        else:
            return False


class HasTau(EventFilter):
    """Keep events with at least one good tau"""

    def passes(self, event):
        return len(event.taus) > 0


############################################################
# AntiVBF filter tools
############################################################

## Overlap removal tools
def getPhotons(event):
    """Returns the photons to be removed from the truth jets"""

    truthPhotons = []

    for p in event.mc:
        if not p.pdgId == 22: continue
        if not p.status == 1: continue
        if not p.fourvect.Pt() >= 15*GeV: continue
        if not abs(p.fourvect.Eta()) <= 2.5: continue
        truthPhotons.append(p.fourvect)

    return truthPhotons


def getElectrons(event):
    """Returns the electrons to be removed from the truth jets"""

    truthElectrons = []

    for p in event.mc:
        if not abs(p.pdgId) == 11: continue
        if not p.status == 1: continue
        if not p.fourvect.Pt() >= 15*GeV: continue
        if not abs(p.fourvect.Eta()) <= 2.5: continue
        truthElectrons.append(p.fourvect)

    return truthElectrons


def getTaus(event):
    """Returns the hadronic taus to be removed from the truth jets"""

    truthTaus = []

    for p in event.mc:
        if not abs(p.pdgId) == 15: continue
        if p.status == 3: continue

        leptonic = False

        tauNeutrino = LorentzVector()
        
        for c in p.iter_children():
            if abs(c.pdgId) == 12: leptonic = True
            if abs(c.pdgId) == 14: leptonic = True
            if abs(c.pdgId) == 15: leptonic = True
            if abs(c.pdgId) == 16: tauNeutrino = c.fourvect

        if leptonic: continue

        visTau = p.fourvect - tauNeutrino

        if not visTau.Pt() >= 15*GeV: continue
        if not abs(visTau.Eta()) <= 2.5: continue

        truthTaus.append(visTau)

    return truthTaus


def truthJetOverlap(truthJet, vectors):
    """ Test the dR overlap of truthJet with all vectors """

    for v in vectors:
        dR = truthJet.DeltaR(v)
        ptRatio = abs(truthJet.Pt() - v.Pt())/v.Pt()
        if dR <= 0.05 and ptRatio <= 0.3:
            return True

    return False

def vbfFilter(event, deta_cut=2.0, mjj_cut=200):
    jets = []

    overlapParticles = getPhotons(event) + getElectrons(event) + getTaus(event)
        
    for jet in event.truthjets:
        if jet.fourvect.Pt() < 15*GeV: continue
        if abs(jet.fourvect.Eta()) > 5.0: continue
        if truthJetOverlap(jet.fourvect, overlapParticles): continue
        jets.append(jet.fourvect)

    passNjets = False

    njets =  len(jets)
    if njets >= 2 : passNjets = True

    passMjj  = False
    passdEta = False

    for i in range(njets):
        for j in range(njets):
            if j > i:
                dEta = abs(jets[i].Eta() - jets[j].Eta())
                Mjj  = (jets[i] + jets[j]).M()
                
                if dEta > deta_cut: passdEta = True
                if Mjj > mjj_cut*GeV: passMjj = True

    if passMjj and passdEta and passNjets:
        return True
    return False


def vbfFilter3(event, deta_cut=2.0, mjj_cut=200):

    overlapParticles = getPhotons(event) + getElectrons(event) + getTaus(event)
        
    for jet in event.truthjets:
        if jet.fourvect.Pt() < 15*GeV: continue
        if abs(jet.fourvect.Eta()) > 5.0: continue
        if truthJetOverlap(jet.fourvect, overlapParticles): continue
        jets.append(jet.fourvect)

    passNjets = False

    njets =  len(jets)
    if njets >= 2 : passNjets = True

    passMjj  = False
    passdEta = False

    for i in event.truthjets:
        for j in event.truthjets:
            for k in event.truthjets:
                if i.Eta() > j.Eta() > k.Eta():
                    Mjj = (i + k).M()
                    dEta = abs(i.Eta() - k.Eta())

                    if dEta > deta_cut: passdEta = True
                    if Mjj > mjj_cut*GeV: passMjj = True

    if passMjj and passdEta and passNjets:
        return True
    return False


class VBFFilter(EventFilter):
    """Keep events that don't pass the VBF filter"""

    def __init__(self, deta_cut, mjj_cut, **kwargs):

        self.deta_cut = deta_cut
        self.mjj_cut
        super(VBFFilter, self).__init__(**kwargs)

    def passes(self, event):
        return vbfFilter(event, self.deta_cut, self.mjj_cut)


class VBFFilter3(EventFilter):
    """Keep events that don't pass the VBF filter"""

    def __init__(self, deta_cut, mjj_cut, **kwargs):

        self.deta_cut = deta_cut
        self.mjj_cut
        super(VBFFilter3, self).__init__(**kwargs)

    def passes(self, event):
        return vbfFilter3(event, self.deta_cut, self.mjj_cut)
            

class AntiVBFFilter(EventFilter):
    """Keep events that don't pass the VBF filter"""

    def passes(self, event):
        return not vbfFilter(event)


                        

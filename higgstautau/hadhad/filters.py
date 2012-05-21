from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
import datetime

"""
See main documentation here:
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauToHH2012Winter
"""

MIN_TAUS = 2


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

        event.taus.select(lambda tau: tau.pt > 20 * GeV)
        return len(event.taus) >= MIN_TAUS


class TauEta(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: abs(tau.eta) < 2.5) # was 2.1
        return len(event.taus) >= MIN_TAUS


class TauJVF(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: True if
                abs(tau.track_eta[tau.leadtrack_idx]) > 2.1 else tau.jet_jvtxf > .5)
        return len(event.taus) >= MIN_TAUS


class Tau1Track3Track(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.numTrack in (1, 3))
        return len(event.taus) >= MIN_TAUS


class TauCharge(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: abs(tau.charge) == 1)
        return len(event.taus) >= MIN_TAUS


class TauIDLoose(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.JetBDTSigLoose)
        return len(event.taus) >= MIN_TAUS


class TauIDMedium(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: tau.JetBDTSigMedium)
        return len(event.taus) >= MIN_TAUS


class TauCrack(EventFilter):

    def passes(self, event):

        event.taus.select(lambda tau: not (1.37 < abs(tau.track_eta[tau.leadtrack_idx]) < 1.52))
        return len(event.taus) >= MIN_TAUS


class TauLArHole(EventFilter):

    def passes(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        event.taus.select(lambda tau:
                not (-0.1 < tau.track_eta[tau.leadtrack_idx] < 1.55
                 and -0.9 < tau.track_phi[tau.leadtrack_idx] < -0.5))
        return len(event.taus) >= MIN_TAUS


class TauLeadSublead(EventFilter):

    def __init__(self, lead = 25*GeV, sublead = 35*GeV, **kwargs):

        super(TauLeadSublead, self).__init__(**kwargs)
        """
        Leading and subleading tau pT thresholds
        """
        self.lead = lead
        self.sublead = sublead

    def passes(self, event):
        # sort in descending order by pT
        event.taus.sort(key=lambda tau: tau.pt, reverse=True)
        # only keep leading two taus
        event.taus.slice(0, 2)
        # Event passes if the highest pT tau is above the leading
        # pT threshold and the next subleading tau pT is above the subleading pT theshold
        return event.taus[0].pt > self.lead and event.taus[1].pt > self.sublead


class TauTriggerMatch(EventFilter):

    def __init__(self,
                 config,
                 datatype,
                 dR=0.2,
                 year=None,
                 **kwargs):

        super(TauTriggerMatch, self).__init__(**kwargs)
        self.config = config
        self.dR = dR
        if year is None:
            year = datetime.datetime.now().year
        year %= 1000

        """
        WARNING: possible bias if matching between MC and data differs
        """
        if datatype == datasets.DATA:
            from ..trigger import utils as triggerutils
            self.triggerutils = triggerutils
            if year == 11:
                self.passes = self.passes_data11
            else:
                raise ValueError("No trigger matching defined for year %d" % year)
        else:
            if year == 11:
                self.passes = self.passes_mc11
            else:
                raise ValueError("No trigger matching defined for year %d" % year)

    def passes_mc11(self, event):

        """
        Matching performed during first skim with CoEPPTrigTool
        """
        event.taus.select(lambda tau: tau.trigger_match_index > -1)
        return len(event.taus) == MIN_TAUS

    def passes_data11(self, event):

        if 177986 <= event.RunNumber <= 187815: # Periods B-K
            trigger = 'EF_tau29_medium1_tau20_medium1'
        elif 188902 <= event.RunNumber <= 191933: # Periods L-M
            trigger = 'EF_tau29T_medium1_tau20T_medium1'
        else:
            raise ValueError("No trigger defined for run %i" % event.RunNumber)

        # erase any previous selection on trigger taus
        # (just to be safe, there should not be any)
        event.taus_EF.reset()
        # get indices of trigger taus associated with this trigger
        trigger_idx = self.triggerutils.get_tau_trigger_obj_idx(
                            self.config,
                            event,
                            trigger)
        # select the trigger taus
        event.taus_EF.select_indices(trigger_idx)
        #print list(event.taus_EF)
        matched_taus = []
        for triggertau in event.taus_EF:
            # find closest matching reco offline tau
            closest_dR = 99999
            closest_tau = None
            closest_idx = -1
            for i, tau in enumerate(event.taus):
                dR = utils.dR(triggertau.eta, triggertau.phi, tau.eta, tau.phi)
                if dR < self.dR and dR < closest_dR:
                    closest_dR = dR
                    closest_tau = tau
                    closest_idx = i
            if closest_tau is not None:
                matched_taus.append(closest_tau)
                # remove match from consideration by future matches (greedy match)
                event.taus.mask_indices([closest_idx])
        # select only the matching offline taus (if any)
        event.taus.reset()
        event.taus.select_indices([tau.index for tau in matched_taus])
        # event passes if at least two taus selected
        return len(event.taus) == MIN_TAUS


class Triggers(EventFilter):
    """
    See lowest unprescaled triggers here:
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled#Taus_electron_muon_MET
    """
    triggers_11 = [
        'EF_tau29_medium1_tau20_medium1',
        'EF_tau29T_medium1_tau20T_medium1'
    ]

    triggers_12 = [
        'EF_tau29Ti_medium1_tau20Ti_medium1',
        'EF_2tau38T_medium1'
    ]

    def __init__(self, datatype, year=None, skim=False, **kwargs):

        if year is None:
            year = datetime.datetime.now().year
        year %= 1000
        if year == 11:
            if datatype == datasets.DATA:
                self.passes = self.passes_data11
            else:
                if skim:
                    self.passes = self.passes_mc11_skim
                else:
                    self.passes = self.passes_mc11
        elif year == 12:
            if datatype == datasets.DATA:
                self.passes = self.passes_data12
            else:
                self.passes = self.passes_mc12
        else:
            raise ValueError("No triggers defined for year %d" % year)
        super(Triggers, self).__init__(**kwargs)

    def passes_mc11(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1_EMULATED
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1_EMULATED
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)

    def passes_mc11_skim(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)

    def passes_data11(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)

    def passes_data12(self, event):
        try:
            return event.EF_tau29Ti_medium1_tau20Ti_medium1 or event.EF_2tau38T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e

    def passes_mc12(self, event):
        try:
            return event.EF_tau29Ti_medium1_tau20Ti_medium1 or event.EF_2tau38T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e


data_triggers = [
    'EF_tau29_medium1_tau20_medium1',
    'EF_tau29T_medium1_tau20T_medium1'
    'EF_tau100_medium',
    'EF_tau125_medium1',
    'EF_xe60_noMu',
    'EF_xe60_tight_noMu',
    'EF_xe60_verytight_noMu',
]

mc_triggers = [
    'EF_tau29_medium1_tau20_medium1',
    'L1_2TAU11_TAU15',
    'EF_tau100_medium',
    'EF_tau125_medium1',
    'EF_xe60_noMu',
    'EF_xe60_tight_noMu',
    'EF_xe60_verytight_noMu'
]


class SkimmingDataTriggers(EventFilter):

    def passes(self, event):
        """
        The lowest un-prescaled trigger can be found in
        https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled

        Periods     Runs            Luminosity      Trigger selection
        B           177986-178109   17.379 pb^-1    tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        D           179710-180481   184.774 pb^-1   tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        E           180614-180776   52.198 pb^-1    tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        F           182013-182519   157.298 pb^-1   tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        G           182726-183462   571.944 pb^-1   tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        H           183544-184169   285.787 pb^-1   tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        I           185353-186493   410.998 pb^-1   tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_noMu
        J           186516-186755   239.713 pb^-1   tau29_medium1_tau20_medium1 OR tau100_medium OR xe60_tight_noMu
        K1-K2       186873-187219                   tau29_medium1_tau20_medium1 OR tau125_medium1 OR xe60_tight_noMu
        K3-K6       187453-187815   683.897 pb^-1   EF_tau29_medium1_tau20_medium1 OR tau29T_medium1_tau20T_medium1 OR tau125_medium1 OR xe60_tight_noMu
        L           188902-190343   1613.055 pb^-1  tau29T_medium1_tau20T_medium1 OR tau125_medium1 OR xe60_verytight_noMu
        M           190608-191933   1172.250 pb^-1  tau29T_medium1_tau20T_medium1 OR tau125_medium1 OR xe60_verytight_noMu
        """
        try:
            if 177986 <= event.RunNumber <= 186493: # Period B-I 1480.36 pb-1
                return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau100_medium or event.EF_xe60_noMu
            elif 186516 <= event.RunNumber <= 186755: # Period J 226.392 pb-1
                return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau100_medium or event.EF_xe60_tight_noMu
            elif 186873 <= event.RunNumber <= 187219: # Period K1-K2 391.09 pb-1
                return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau125_medium1 or event.EF_xe60_tight_noMu
            elif 187453 <= event.RunNumber <= 187815: # Period K3-K6 170.616 pb-1
                return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau29T_medium1_tau20T_medium1 or event.EF_tau125_medium1 or event.EF_xe60_tight_noMu
            elif 188902 <= event.RunNumber <= 191933: # Period L-M 2392.85 pb-1
                return event.EF_tau29T_medium1_tau20T_medium1 or event.EF_tau125_medium1 or event.EF_xe60_verytight_noMu
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)


class SkimmingMCTriggers(EventFilter):

    def passes(self, event):
        """
        OR of all triggers above

        EF_tau29T_medium1_tau20T_medium1 is not available in new MC11 so we use the equivalent:
        EF_tau29_medium1_tau20_medium1 && L1_2TAU11_TAU15
        """
        return event.EF_tau29_medium1_tau20_medium1 or event.EF_tau100_medium or event.EF_xe60_noMu or \
               event.EF_xe60_tight_noMu or event.EF_tau125_medium1 or event.EF_xe60_verytight_noMu


class JetCrackVeto(EventFilter):

    def passes(self, event):

        for jet in event.jets:
            if jet.pt <= 20 * GeV: continue
            if 1.3 < abs(jet.eta) < 1.7: return False
        return True


class ElectronVeto(EventFilter):

    def passes(self, event):

       for el in event.electrons:
           pt = el.cl_E / cosh(el.tracketa)
           if pt <= 15 * GeV: continue
           if not ((abs(el.tracketa) < 1.37) or (1.52 < abs(el.tracketa) < 2.47)): continue
           if el.author not in (1, 3): continue
           if not abs(el.charge) == 1: continue
           if el.mediumPP != 1: continue
           if (el.OQ & 1446) != 0: continue
           return False
       return True


def muon_has_good_track(muon):

    blayer = (muon.expectBLayerHit == 0) or (muon.nBLHits > 0)
    pix = muon.nPixHits + muon.nPixelDeadSensors > 1
    sct = muon.nSCTHits + muon.nSCTDeadSensors >= 6
    holes = muon.nPixHoles + muon.nSCTHoles < 3
    trt = False

    if abs(muon.eta) < 1.9:
        trt = (muon.nTRTHits + muon.nTRTOutliers > 5) and \
              (muon.nTRTOutliers < 0.9 * (muon.nTRTHits + muon.nTRTOutliers))
    else:
        trt = (muon.nTRTHits + muon.nTRTOutliers <= 5) or \
              (muon.nTRTOutliers < 0.9 * (muon.nTRTHits + muon.nTRTOutliers))

    return blayer and pix and sct and holes and trt


class MuonVeto(EventFilter):

    def passes(self, event):

       for muon in event.muons:
           if muon.pt <= 10 * GeV: continue
           if abs(muon.eta) >= 2.5: continue
           if muon.loose != 1: continue
           if not muon_has_good_track(muon): continue
           return False

       return True

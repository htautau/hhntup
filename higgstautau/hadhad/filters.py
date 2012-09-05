from rootpy.tree.filtering import *
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *


class TauLeadSublead(EventFilter):

    def __init__(self, lead=35*GeV, sublead=25*GeV, **kwargs):

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

    # put triggers in descending order of pT thresholds
    # put thresholds in descending order
    triggers_11 = [
        ('EF_tau29_medium1_tau20_medium1', (29, 20)),
        ('EF_tau29T_medium1_tau20T_medium1', (29, 20))
    ]

    triggers_12 = [
        ('EF_2tau38T_medium1', (38,)),
        ('EF_tau29Ti_medium1_tau20Ti_medium1', (29, 20))
    ]

    def __init__(self,
                 config,
                 datatype,
                 year,
                 dR=0.2,
                 skim=False,
                 tree=None,
                 num_taus=2,
                 min_taus=None,
                 **kwargs):

        super(TauTriggerMatch, self).__init__(**kwargs)
        self.config = config
        self.dR = dR
        year %= 1000
        self.skim = skim
        self.tree = tree
        self.num_taus = num_taus
        self.min_taus = min_taus

        """
        WARNING: possible bias if matching between MC and data differs
        """
        if datatype == datasets.DATA:
            from ..trigger import utils as triggerutils
            self.triggerutils = triggerutils
            if year == 11:
                self.passes = self.passes_data11
            elif year == 12:
                self.passes = self.passes_data12
            else:
                raise ValueError("No data trigger matching defined for year %d" % year)
        else: # MC
            if year == 11:
                self.passes = self.passes_mc11
            elif year == 12:
                from ..trigger import utils as triggerutils
                self.triggerutils = triggerutils
                self.passes = self.passes_mc12
            else:
                raise ValueError("No MC trigger matching defined for year %d" % year)

    def passes_mc11(self, event):

        """
        Matching performed during first skim with CoEPPTrigTool
        """
        event.taus.select(lambda tau: tau.trigger_match_index > -1)
        if self.min_taus is not None:
            return len(event.taus) >= self.min_taus
        else:
            return len(event.taus) == self.num_taus

    def passes_mc12(self, event):

        self.match(event, self.triggers_12)
        if self.min_taus is not None:
            return len(event.taus) >= self.min_taus
        else:
            return len(event.taus) == self.num_taus

    def passes_data11(self, event):

        if 177986 <= event.RunNumber <= 187815: # Periods B-K
            trigger = 'EF_tau29_medium1_tau20_medium1'
        elif 188902 <= event.RunNumber <= 191933: # Periods L-M
            trigger = 'EF_tau29T_medium1_tau20T_medium1'
        else:
            raise ValueError("No trigger defined for run %i" % event.RunNumber)
        # TODO: clean up trigger config (no hardcoded values...)
        self.match(event, [(trigger, (29, 20))])
        if self.min_taus is not None:
            return len(event.taus) >= self.min_taus
        else:
            return len(event.taus) == self.num_taus

    def passes_data12(self, event):

        self.match(event, self.triggers_12)
        if self.min_taus is not None:
            return len(event.taus) >= self.min_taus
        else:
            return len(event.taus) == self.num_taus

    def match(self, event, triggers):
        """
        triggers: list of triggers (and thresholds) that are ORed.
        Order them by priority i.e. put triggers with higher pT
        thesholds before triggers with lower thresholds.
        """
        matched_taus = []
        matches = {}
        for trigger, thresholds in triggers:
            # get indices of trigger taus associated with this trigger
            trigger_idx = self.triggerutils.get_tau_trigger_obj_idx(
                self.config,
                event,
                trigger)
            # erase any previous selection on trigger taus
            # (just to be safe, there should not be any)
            event.taus_EF.reset()
            # select the trigger taus
            event.taus_EF.select_indices(trigger_idx)
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
                    if self.skim:
                        for threshold in thresholds:
                            if triggertau.pt > threshold * GeV:
                                thresh = threshold
                                break
                        matches[closest_tau.index] = (triggertau.index, thresh)
                    matched_taus.append(closest_tau)
                    # remove match from consideration by future matches (greedy match)
                    event.taus.mask_indices([closest_idx])
        # select only the matching offline taus (if any)
        event.taus.reset()
        if self.skim:
            self.tree.tau_trigger_match_index.clear()
            self.tree.tau_trigger_match_thresh.clear()
            for i in xrange(event.tau_n):
                if i in matches:
                    idx, thresh = matches[i]
                    self.tree.tau_trigger_match_index.push_back(idx)
                    self.tree.tau_trigger_match_thresh.push_back(thresh)
                else:
                    self.tree.tau_trigger_match_index.push_back(-1)
                    self.tree.tau_trigger_match_thresh.push_back(0)
        event.taus.select_indices([tau.index for tau in matched_taus])


class Triggers(EventFilter):
    """
    See lowest unprescaled triggers here:
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled#Taus_electron_muon_MET
    """
    def __init__(self, year, **kwargs):

        if year == 2011:
            self.passes = self.passes_11
        elif year == 2012:
            self.passes = self.passes_12
        else:
            raise ValueError("No triggers defined for year %d" % year)
        super(Triggers, self).__init__(**kwargs)

    def passes_11(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" % event.RunNumber)

    def passes_12(self, event):
        try:
            return event.EF_tau29Ti_medium1_tau20Ti_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (event.RunNumber, e)
            raise e


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


from ..filters import muon_has_good_track

class MuonVeto(EventFilter):

    def __init__(self, year, **kwargs):

        self.year = year
        super(MuonVeto, self).__init__(**kwargs)

    def passes(self, event):

       for muon in event.muons:
           if muon.pt <= 10 * GeV:
               continue
           if abs(muon.eta) >= 2.5:
               continue
           if muon.loose != 1:
               continue
           if not muon_has_good_track(muon, self.year):
               continue
           return False
       return True

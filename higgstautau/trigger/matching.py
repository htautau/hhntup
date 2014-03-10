from rootpy.tree.filtering import EventFilter

from math import *

from . import utils as triggerutils
from . import log; log = log[__name__]
from .. import utils
from ..units import GeV
from .. import datasets


class TauTriggerMatchIndex(EventFilter):
    """
    Match reco taus to trigger taus. If there are more than two EF taus, take
    the leading two.
    """
    def __init__(self,
                 config,
                 datatype,
                 year,
                 dR=0.2,
                 passthrough=False,
                 **kwargs):

        super(TauTriggerMatchIndex, self).__init__(
                passthrough=passthrough,
                **kwargs)

        if not passthrough:
            self.config = config
            self.dR = dR
            year %= 1000
            self.year = year
            """
            WARNING: possible bias if matching between MC and data differs
            """
            if datatype == datasets.DATA:
                if year == 11:
                    self.passes = self.passes_data11
                elif year == 12:
                    self.passes = self.passes_data12
                else:
                    raise ValueError(
                            "No data trigger matching defined for year %d" %
                            year)
            elif datatype == datasets.MC:
                if year == 11:
                    self.passes = self.passes_mc11
                elif year == 12:
                    self.passes = self.passes_mc12
                else:
                    raise ValueError(
                            "No MC trigger matching defined for year %d" %
                            year)
            else:
                raise ValueError(
                        "No trigger matching defined for datatype %d" %
                        datatype)

    def passes_mc11(self, event):
        """
        Matching performed during trigger emulation with CoEPPTrigTool
        """
        event.taus.select(lambda tau: tau.trigger_match_index > -1)
        self.match_sanity(event)
        return len(event.taus) >= 2

    def passes_data11(self, event):

        if 177986 <= event.RunNumber <= 187815: # Periods B-K
            trigger = 'EF_tau29_medium1_tau20_medium1'
        elif 188902 <= event.RunNumber <= 191933: # Periods L-M
            trigger = 'EF_tau29T_medium1_tau20T_medium1'
        else:
            raise ValueError("No trigger defined for run %i" % event.RunNumber)
        self.match_index(event, trigger)
        event.taus.select(lambda tau: tau.trigger_match_index > -1)
        self.match_sanity(event)
        return len(event.taus) >= 2

    def passes_mc12(self, event):

        self.match_index(event, 'EF_tau29Ti_medium1_tau20Ti_medium1')
        event.taus.select(lambda tau: tau.trigger_match_index > -1)
        self.match_sanity(event)
        return len(event.taus) >= 2

    def passes_data12(self, event):

        self.match_index(event, 'EF_tau29Ti_medium1_tau20Ti_medium1')
        event.taus.select(lambda tau: tau.trigger_match_index > -1)
        self.match_sanity(event)
        return len(event.taus) >= 2

    def match_index(self, event, trigger):

        # get indices of trigger taus associated with this trigger
        trigger_idx = triggerutils.get_tau_trigger_obj_idx(
            self.config,
            event,
            trigger)
        # trigger_idx can contain 3 indices
        # will need to take the leading two

        taus = list(event.taus)

        # for each EF tau find closest matching reco tau
        for EF_idx in trigger_idx:
            trigger_tau = event.taus_EF.getitem(EF_idx)
            closest_dR = 99999
            closest_tau = None
            for tau in taus:
                dR = utils.dR(
                        trigger_tau.eta, trigger_tau.phi,
                        tau.eta, tau.phi)
                if dR < self.dR and dR < closest_dR:
                    closest_dR = dR
                    closest_tau = tau
            if closest_tau is not None:
                closest_tau.trigger_match_index = EF_idx
                # remove match from future matches (greedy match)
                taus.remove(closest_tau)

    def match_sanity(self, event):

        # sanity check
        if len(event.taus) < 3:
            return
        print '-' * 20
        print "Run: %d, Event %d" % (event.RunNumber, event.EventNumber)
        fmt = "%d:   match index: %d    reco pT: %f    EF pT: %f"
        for i, tau in enumerate(event.taus):
            print fmt % (i, tau.trigger_match_index, tau.pt,
                    event.taus_EF.getitem(tau.trigger_match_index).pt)
            print "dR with all other taus:"
            for j, tau2 in enumerate(event.taus):
                if tau2 != tau:
                    print j, utils.dR(
                            tau.eta, tau.phi,
                            tau2.eta, tau2.phi)
            print '='


class TauTriggerMatchThreshold(EventFilter):
    """
    Match previously matched reco taus to thresholds of the trigger
    """
    def __init__(self,
                 datatype,
                 tree,
                 passthrough=False,
                 **kwargs):

        self.datatype = datatype
        self.tree = tree
        super(TauTriggerMatchThreshold, self).__init__(
                passthrough=passthrough,
                **kwargs)

    def passes(self, event):

        if self.datatype == datasets.EMBED:
            assert len(event.taus) == 2
            assert event.taus[0].pt >= event.taus[1].pt
            # taus are already sorted in descending order by pT by TauLeadSublead
            tau1, tau2 = event.taus
            tau1.trigger_match_thresh = 29
            tau2.trigger_match_thresh = 20
        else:
            self.match_threshold(event, (29, 20))
        return True

    def match_threshold(self, event, thresholds):
        """
        thresholds must be in descending order

        TODO: Use the info stored in the D3PD for 2012:
        trig_EF_tau_EF_tau29Ti_medium1_tau20Ti_medium1

        if(trig_EF_tau_EF_tau20Ti_medium1==1) tau_trigger_match_thresh = 20
        if(trig_EF_tau_EF_tau29Ti_medium1==1) tau_trigger_match_thresh = 29

        """
        assert len(event.taus) == 2
        assert len(thresholds) == 2

        # assume only matched taus remain in event.taus
        taus = [(tau, event.taus_EF.getitem(tau.trigger_match_index)) for
                tau in event.taus]

        # sort by pT of EF tau
        taus = sorted(taus, key=lambda tau: tau[1].pt, reverse=True)

        # assign thresholds in descending order
        for i in xrange(len(taus)):
            taus[i][0].trigger_match_thresh = thresholds[i]
            # sanity check THIS SOMETIMES FAILS!
            if taus[i][1].pt < thresholds[i] * GeV:
                log.warning("EF pT %f less than trigger threshold %f" % (
                    taus[i][1].pt, thresholds[i] * GeV))
                self.tree.tau_trigger_match_error = True

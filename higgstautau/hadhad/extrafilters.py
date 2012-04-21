from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *

from ..trigger import utils as triggerutils

MIN_TAUS = 2


class TauTriggerMatch(EventFilter):

    def __init__(self, config, dR=0.2, **kwargs):

        super(TauTriggerMatch, self).__init__(**kwargs)
        self.config = config
        self.dR = dR

    def passes(self, event):

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
        trigger_idx = triggerutils.get_tau_trigger_obj_idx(
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



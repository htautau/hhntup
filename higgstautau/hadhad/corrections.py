import os
import ROOT
from rootpy.tree.filtering import *
from externaltools import TauTriggerCorrections


class TauTriggerEfficiency(EventFilter):
    """
    This filter should only be run on MC
    currently only supports 2011
    """

    def __init__(self, year, **kwargs):

        self.year = year
        base = TauTriggerCorrections.RESOURCE_PATH
        self.correct_20 = ROOT.TauTriggerCorrections(os.path.join(base,
            'triggerSF_EF_tau20_medium1.root'))
        self.correct_29 = ROOT.TauTriggerCorrections(os.path.join(base,
            'triggerSF_EF_tau29_medium1.root'))
        self.correct_20T = ROOT.TauTriggerCorrections(os.path.join(base,
            'triggerSF_EF_tau20T_medium1.root'))
        self.correct_29T = ROOT.TauTriggerCorrections(os.path.join(base,
            'triggerSF_EF_tau29T_medium1.root'))
        super(TauTriggerEfficiency, self).__init__(**kwargs)

    def passes(self, event):

        thresh = []
        for tau in event.taus:
            if tau.trigger_match_thresh == 20:
                if event.RunNumber >= 188902:
                     corr = self.correct_20T
                else:
                     corr = self.correct_20
                thresh.append(20)
            elif tau.trigger_match_thresh == 29:
                if event.RunNumber >= 188902:
                    corr = self.correct_29T
                else:
                    corr = self.correct_29
                thresh.append(29)
            else:
                raise ValueError("trigger match thresh of %d is not understood"
                        % tau.trigger_match_thresh)

            tau.trigger_scale_factor = corr.getSF(tau.pt, 0)
            tau.trigger_scale_factor_high = corr.getSF(tau.pt, 1)
            tau.trigger_scale_factor_low = corr.getSF(tau.pt, -1)

        if len(set(thresh)) != 2 or len(thresh) != 2:
            raise Exception("there must be exactly two unique trigger match"
                    " thresholds (29, 20). Got: %s" % str(thresh))
        return True

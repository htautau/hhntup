import os
from operator import itemgetter
import ROOT
from rootpy.tree.filtering import *
from atlastools.units import GeV
from externaltools import TauTriggerCorrections
from atlastools import datasets


class TauTriggerEfficiency(EventFilter):
    """
    This filter should only be run on MC
    currently only supports 2011
    """

    def __init__(self, year, datatype, **kwargs):

        self.year = year % 1000
        self.datatype = datatype
        base = TauTriggerCorrections.RESOURCE_PATH

        if self.year == 11:
            if datatype == datasets.MC:
                self.correct_20 = ROOT.TauTriggerCorrections(os.path.join(base,
                    'triggerSF_EF_tau20_medium1.root'))
                self.correct_29 = ROOT.TauTriggerCorrections(os.path.join(base,
                    'triggerSF_EF_tau29_medium1.root'))
                self.correct_20T = ROOT.TauTriggerCorrections(os.path.join(base,
                    'triggerSF_EF_tau20T_medium1.root'))
                self.correct_29T = ROOT.TauTriggerCorrections(os.path.join(base,
                    'triggerSF_EF_tau29T_medium1.root'))
            elif datatype == datasets.EMBED:
                # use the 3D param
                self.correct_20 = ROOT.TauTriggerCorrections()
                self.correct_20.loadInputFile(
                        "../root/triggerSF_wmcpara_EF_tau20_medium1.root",
                        "1P3P","BDTm" )
                # separate for medium and tight
            else:
                raise ValueError(
                    "No trigger efficiency corrections defined for datatype %d"
                    % datatype)
            self.passes = self.passes_11
        else:
            raise ValueError(
                "No trigger efficiency corrections defined for year %d" % year)
        super(TauTriggerEfficiency, self).__init__(**kwargs)

    def passes_11(self, event):

        assert len(event.taus) == 2

        idx = [tau.trigger_match_index for tau in event.taus]
        assert len(set(idx)) == 2

        # fix trigger threshold association
        taus = [(tau, event.taus_EF.getitem(tau.trigger_match_index)) for
                tau in event.taus]
        # sort by pT of EF tau
        taus = sorted(taus, key=lambda tau: tau[1].pt, reverse=True)

        taus[0][0].trigger_match_thresh = 29
        taus[1][0].trigger_match_thresh = 20

        assert taus[0][1].pt > 29 * GeV
        assert taus[1][1].pt > 20 * GeV

        if event.RunNumber >= 188902:
            assert taus[0][1].EF_tau29T_medium1 == 1
            assert taus[1][1].EF_tau20T_medium1 == 1
        else:
            assert taus[0][1].EF_tau29_medium1 == 1
            assert taus[1][1].EF_tau20_medium1 == 1

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

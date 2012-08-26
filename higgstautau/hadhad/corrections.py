import os
from operator import itemgetter
import ROOT
from rootpy.tree.filtering import *
from atlastools.units import GeV
from externaltools import TauTriggerCorrections
from atlastools import datasets

from ..filters import pileup_vertex_selection


class TauTriggerEfficiency(EventFilter):
    """
    This filter should only be run on MC
    currently only supports 2011
    """

    def __init__(self, year, datatype, tes_systematic=False, **kwargs):

        self.year = year % 1000
        self.datatype = datatype
        self.tes_systematic = tes_systematic
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
                self.passes = self.passes_11_mc
            elif datatype == datasets.EMBED:
                # use the 3D param
                self.corrections = {}
                for thresh in (29, 20):
                    self.corrections[thresh] = {}
                    for tight, bdt_label  in zip((0, 1), ('BDTm', 'BDTt')):
                        self.corrections[thresh][tight] = {}
                        for period, is_late in zip(('', 'T'), (False, True)):
                            tool = ROOT.TauTriggerCorrections()
                            tool.loadInputFile(
                                os.path.join(base,
                                "triggerSF_wmcpara_EF_tau%d%s_medium1.root" % (
                                    thresh, period)),
                                "1P3P", bdt_label)
                            self.corrections[thresh][tight][is_late] = tool
                self.passes = self.passes_11_embed
            else:
                raise ValueError(
                    "No trigger efficiency corrections defined for datatype %d"
                    % datatype)
        else:
            raise ValueError(
                "No trigger efficiency corrections defined for year %d" % year)
        super(TauTriggerEfficiency, self).__init__(**kwargs)

    def passes_11_mc(self, event):

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

            # correct for TES variations (only on nominal)
            # pt_nominal should equal pt when TES is not applied
            if self.tes_systematic:
                print "%f %f" % (tau.pt, tau.pt_nominal)
                tau.trigger_scale_factor = (corr.getSF(tau.pt, 0) *
                    corr.getEffMC(tau.pt, 0) / corr.getEffMC(tau.pt_nominal, 0))
            else:
                tau.trigger_scale_factor = corr.getSF(tau.pt, 0)
            tau.trigger_scale_factor_high = corr.getSF(tau.pt, 1)
            tau.trigger_scale_factor_low = corr.getSF(tau.pt, -1)

        if len(set(thresh)) != 2 or len(thresh) != 2:
            raise Exception("there must be exactly two unique trigger match"
                    " thresholds (29, 20). Got: %s" % str(thresh))
        return True

    def passes_11_embed(self, event):

        assert len(event.taus) == 2
        # taus are already sorted in descending order by pT by TauLeadSublead

        npileup_vtx = len([vtx for vtx in event.vertices
                           if pileup_vertex_selection(vtx)])

        for tau, thresh in zip(event.taus, (29, 20)):

            corr = self.corrections[thresh][tau.JetBDTSigTight][event.RunNumber >= 188902]
            tau.trigger_scale_factor = corr.get3DMCEff(
                        tau.pt, tau.eta,
                        npileup_vtx, 0)
            tau.trigger_scale_factor_high = corr.get3DMCEff(
                        tau.pt, tau.eta,
                        npileup_vtx, 1)
            tau.trigger_scale_factor_low = corr.get3DMCEff(
                        tau.pt, tau.eta,
                        npileup_vtx, -1)

        return True

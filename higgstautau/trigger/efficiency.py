import os
import math
from operator import itemgetter
import ROOT
from rootpy.tree.filtering import *
from atlastools.units import GeV
from atlastools import datasets

from ..filters import pileup_vertex_selection
from ..tauid import IDLOOSE, IDMEDIUM, IDTIGHT
from . import log; log = log[__name__]


class TauTriggerEfficiency(EventFilter):
    """
    Determine the trigger efficiency corrections
    """
    @classmethod
    def get_period_2012(cls, run):

        if run <= 201556:
            return 'periodA'
        elif run <= 209025:
            return 'periodBD'
        elif run <= 216432:
            return 'periodEM'
        # FAIL
        return None

    def __init__(self,
            year,
            datatype,
            tree,
            tes_systematic=False,
            passthrough=False,
            **kwargs):

        if not passthrough:
            self.year = year % 1000
            self.datatype = datatype
            self.tes_systematic = tes_systematic
            self.tree = tree

            if self.year == 11:
                from externaltools.bundle_2011 import TauTriggerCorrections
                base = TauTriggerCorrections.RESOURCE_PATH

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
                        for wplevel in ('BDTl', 'BDTm', 'BDTtEVm'):
                            self.corrections[thresh][wplevel] = {}
                            for prong in (1, 3):
                                self.corrections[thresh][wplevel][prong] = {}
                                for period, is_late in zip(('', 'T'), (False, True)):
                                    tool = ROOT.TauTriggerCorrections()
                                    tool.loadInputFile(
                                        os.path.join(base,
                                        "triggerSF_wmcpara_EF_tau%d%s_medium1.root" % (
                                            thresh, period)),
                                        "%dP" % prong, wplevel)
                                    self.corrections[thresh][wplevel][prong][is_late] = tool
                    self.passes = self.passes_11_embed
                else:
                    raise ValueError("datatype is not EMBED or MC")

            elif self.year == 12:
                from externaltools.bundle_2012 import TrigTauEfficiency
                base = os.path.join(TrigTauEfficiency.RESOURCE_PATH,
                        'benchmark_menu')

                if datatype in (datasets.MC, datasets.EMBED):
                    self.ttc_20 = ROOT.TrigTauEfficiency()
                    status = self.ttc_20.loadInputFile(os.path.join(base,
                        'triggerSF_EF_tau20Ti_medium1.root'))
                    if status != 0:
                        raise RuntimeError(
                            'could not load triggerSF_EF_tau20Ti_medium1.root')
                    self.ttc_29 = ROOT.TrigTauEfficiency()
                    status = self.ttc_29.loadInputFile(os.path.join(base,
                        'triggerSF_EF_tau29Ti_medium1.root'))
                    if status != 0:
                        raise RuntimeError(
                            'could not load triggerSF_EF_tau29Ti_medium1.root')
                    if datatype == datasets.MC:
                        self.passes = self.passes_12_mc
                    else:
                        self.passes = self.passes_12_embed
                else:
                    raise ValueError("datatype is not EMBED or MC")
            else:
                raise ValueError(
                    "No trigger efficiency corrections defined for year %d" % year)

            if datatype == datasets.MC:
                log.info("will apply MC trigger scale factors")
            else:
                log.info("will apply data trigger efficiencies")

        super(TauTriggerEfficiency, self).__init__(
                passthrough=passthrough,
                **kwargs)

    def passes_11_mc(self, event):

        """
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
        """

        """
        if self.tree.RunNumber >= 188902:
            assert taus[0][1].EF_tau29T_medium1 == 1
            assert taus[1][1].EF_tau20T_medium1 == 1
        else:
            assert taus[0][1].EF_tau29_medium1 == 1
            assert taus[1][1].EF_tau20_medium1 == 1
        """

        thresh = []
        for tau in event.taus:

            # only correct matched taus
            if not tau.matched:
                continue

            if tau.trigger_match_thresh == 20:
                if self.tree.RunNumber >= 188902:
                     corr = self.correct_20T
                else:
                     corr = self.correct_20
                thresh.append(20)
            elif tau.trigger_match_thresh == 29:
                if self.tree.RunNumber >= 188902:
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
                #print "%f %f" % (tau.pt, tau.pt_nominal)
                tau.trigger_scale_factor = abs(
                    corr.getSF(tau.pt, 0) *
                    corr.getMCEff(tau.pt, 0) /
                    corr.getMCEff(tau.pt_nominal, 0))
            else:
                tau.trigger_scale_factor = abs(corr.getSF(tau.pt, 0))
            tau.trigger_scale_factor_high = abs(corr.getSF(tau.pt, 1))
            tau.trigger_scale_factor_low = abs(corr.getSF(tau.pt, -1))

        return True

    def get_id_11_embed(self, tau):

        if tau.id == IDLOOSE:
            return 'BDTl'
        elif tau.id == IDMEDIUM:
            return 'BDTm'
        elif tau.id == IDTIGHT:
            return 'BDTtEVm'
        raise ValueError("tau is not loose, medium, or tight")

    def passes_11_embed(self, event):

        assert len(event.taus) == 2
        assert event.taus[0].pt > event.taus[1].pt
        # taus are already sorted in descending order by pT by TauLeadSublead

        npileup_vtx = len([vtx for vtx in event.vertices
                           if pileup_vertex_selection(vtx)])

        for tau, thresh in zip(event.taus, (29, 20)):

            # only correct matched taus
            if not tau.matched:
                continue

            wpflag = self.get_id_11_embed(tau)

            if tau.numTrack > 1:
                prong = 3
            else:
                prong = 1

            corr = self.corrections[thresh][wpflag][prong][self.tree.RunNumber >= 188902]
            tau.trigger_scale_factor = abs(corr.get3DMCEff(
                tau.pt, tau.eta,
                npileup_vtx, 0))
            tau.trigger_scale_factor_high = abs(corr.get3DMCEff(
                tau.pt, tau.eta,
                npileup_vtx, 1))
            tau.trigger_scale_factor_low = abs(corr.get3DMCEff(
                tau.pt, tau.eta,
                npileup_vtx, -1))

        return True

    def get_id_12(self, tau):

        if tau.id == IDLOOSE:
            return 'BDTl'
        elif tau.id == IDMEDIUM:
            return 'BDTm'
        elif tau.id == IDTIGHT:
            return 'BDTt'
        raise ValueError("tau is not loose, medium, or tight")

    def passes_12_mc(self, event):

        period = self.get_period_2012(self.tree.RunNumber)
        thresh = []

        # new tool only accepts EVnone
        eveto = 'EVnone'

        for tau in event.taus:

            # only correct matched taus
            if not tau.matched:
                continue

            wpflag = self.get_id_12(tau)

            if tau.numTrack > 1:
                prong = '3p'
            else:
                prong = '1p'

            if tau.trigger_match_thresh == 20:
                ttc = self.ttc_20
                thresh.append(20)
            elif tau.trigger_match_thresh == 29:
                ttc = self.ttc_29
                if tau.pt < 35 * GeV:
                    log.warning("error in trigger scale factors")
                    # should be rare (swapping of reco vs trigger pT)
                    self.tree.tau_trigger_match_error = True
                thresh.append(29)
            else:
                raise ValueError("trigger match thresh of %d is not understood"
                        % tau.trigger_match_thresh)

            if self.tes_systematic:
                # correct for TES variations (only on nominal)
                # pt_nominal should equal pt when TES is not applied
                #print "%f %f" % (tau.pt, tau.pt_nominal)
                try:
                    sf = abs(ttc.getSF(tau.pt, tau.eta, 0, period, prong, wpflag, eveto) *
                        ttc.getMCEff(tau.pt, tau.eta, 0, period, prong, wpflag, eveto) /
                        ttc.getMCEff(tau.pt_nominal, tau.eta, 0, period, prong, wpflag, eveto))
                except ZeroDivisionError:
                    log.warning(
                        "division by zero in trigger scale factors: using 1.")
                    sf = 1.
                if math.isinf(sf) or math.isnan(sf):
                    log.warning("trigger data efficiency is infinite or NaN! Using 0.")
                    tau.trigger_scale_factor = 0.
                    continue
                tau.trigger_scale_factor = sf

            else:
                sf = abs(ttc.getSF(tau.pt, tau.eta, 0, period, prong, wpflag, eveto))
                if math.isinf(sf) or math.isnan(sf):
                    log.warning(
                        "trigger efficiency SF is infinite or NaN! Using 0. "
                        "pt: %f, eta: %f, mode: 0, period: %s, prong: %s, wp: %s, eveto: %s"
                        " trigger: %d"
                        % (tau.pt, tau.eta,
                           period, prong,
                           wpflag, eveto,
                           tau.trigger_match_thresh))
                    tau.trigger_scale_factor = 0.
                    tau.trigger_scale_factor_high = 0.
                    tau.trigger_scale_factor_low = 0.
                    continue

                tau.trigger_scale_factor = sf

                # MC stat uncert
                mc_stat_up = abs(ttc.getSF(tau.pt, tau.eta, 1, period, prong, wpflag, eveto))
                mc_stat_dn = abs(ttc.getSF(tau.pt, tau.eta, -1, period, prong, wpflag, eveto))
                # Data stat uncert
                data_stat_up = abs(ttc.getSF(tau.pt, tau.eta, 2, period, prong, wpflag, eveto))
                data_stat_dn = abs(ttc.getSF(tau.pt, tau.eta, -2, period, prong, wpflag, eveto))
                # Systematic uncert
                sys_up = abs(ttc.getSF(tau.pt, tau.eta, 3, period, prong, wpflag, eveto))
                sys_dn = abs(ttc.getSF(tau.pt, tau.eta, -3, period, prong, wpflag, eveto))

                sigma_up = math.sqrt(math.pow(mc_stat_up, 2.) +
                                     math.pow(data_stat_up, 2.) +
                                     math.pow(sys_up, 2.))
                sigma_dn = math.sqrt(math.pow(mc_stat_dn, 2.) +
                                     math.pow(data_stat_dn, 2.) +
                                     math.pow(sys_dn, 2.))

                tau.trigger_scale_factor_high = sf + sigma_up
                tau.trigger_scale_factor_low = sf - sigma_dn

        return True

    def passes_12_embed(self, event):

        assert len(event.taus) == 2
        assert event.taus[0].pt >= event.taus[1].pt
        # taus are already sorted in descending order by pT by TauLeadSublead

        # new tool only accepts EVnone
        eveto = 'EVnone'

        for tau, ttc in zip(event.taus, (self.ttc_29, self.ttc_20)):

            # only correct matched taus
            if not tau.matched:
                continue

            if tau.numTrack > 1:
                prong = '3p'
            else:
                prong = '1p'

            wpflag = self.get_id_12(tau)

            sf = abs(ttc.getDataEff(tau.pt, tau.eta, 0, self.tree.RunNumber, prong, wpflag, eveto))

            #if self.tree.RunNumber == 207528 and event.EventNumber == 2183594:
            #    print sf, tau.pt, tau.eta, prong, wpflag, eveto

            if math.isinf(sf) or math.isnan(sf):
                log.warning("trigger data efficiency is infinite or NaN! Using 0.")
                tau.trigger_scale_factor = 0.
                tau.trigger_scale_factor_high = 0.
                tau.trigger_scale_factor_low = 0.
                continue

            tau.trigger_scale_factor = sf

            # Data stat uncert
            data_stat_up = abs(ttc.getDataEff(tau.pt, tau.eta, 1, self.tree.RunNumber, prong, wpflag, eveto))
            data_stat_dn = abs(ttc.getDataEff(tau.pt, tau.eta, -1, self.tree.RunNumber, prong, wpflag, eveto))
            # Systematic uncert
            sys_up = abs(ttc.getDataEff(tau.pt, tau.eta, 2, self.tree.RunNumber, prong, wpflag, eveto))
            sys_dn = abs(ttc.getDataEff(tau.pt, tau.eta, -2, self.tree.RunNumber, prong, wpflag, eveto))

            sigma_up = math.sqrt(math.pow(data_stat_up, 2.) +
                                 math.pow(sys_up, 2.))
            sigma_dn = math.sqrt(math.pow(data_stat_dn, 2.) +
                                 math.pow(sys_dn, 2.))

            tau.trigger_scale_factor_high = sf + sigma_up
            tau.trigger_scale_factor_low = sf - sigma_dn

        return True

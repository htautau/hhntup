import os
import math
from math import sqrt
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
        assert len(event.taus) == 2
        assert event.taus[0].pt >= event.taus[1].pt
        # taus are already sorted in descending order by pT by TauLeadSublead

        if self.tree.RunNumber >= 188902:
            tools = [self.correct_29T, self.correct_20T]
        else:
            tools = [self.correct_29, self.correct_20]

        for tau, corr in zip(event.taus, tools):

            # only correct matched taus
            if not tau.matched:
                continue

            # correct for TES variations (only on nominal)
            # pt_nominal should equal pt when TES is not applied
            if self.tes_systematic:
                #print "%f %f" % (tau.pt, tau.pt_nominal)
                sf = abs(
                    corr.getSF(tau.pt, 0) *
                    corr.getMCEff(tau.pt, 0) /
                    corr.getMCEff(tau.pt_nominal, 0))
            else:
                sf = abs(corr.getSF(tau.pt, 0))
            tau.trigger_sf = sf
            tau.trigger_sf_high = abs(corr.getSF(tau.pt, 1))
            tau.trigger_sf_low = abs(corr.getSF(tau.pt, -1))
            tau.trigger_sf_mc_stat_high = sf
            tau.trigger_sf_mc_stat_low = sf
            tau.trigger_sf_data_stat_high = sf
            tau.trigger_sf_data_stat_low = sf
            tau.trigger_sf_sys_high = sf
            tau.trigger_sf_sys_low = sf

            eff = corr.getDataEff(tau.pt, 0)
            eff_errup = corr.getDataEff(tau.pt, 1)
            eff_errdn = corr.getDataEff(tau.pt, -1)

            tau.trigger_eff = eff
            tau.trigger_eff_high = eff + eff_errup
            tau.trigger_eff_low = eff - eff_errdn
            tau.trigger_eff_stat_high = eff
            tau.trigger_eff_stat_low = eff
            tau.trigger_eff_sys_high = eff
            tau.trigger_eff_sys_low = eff

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
        assert event.taus[0].pt >= event.taus[1].pt
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

            # TODO: account for TES systematic?
            corr = self.corrections[thresh][wpflag][prong][self.tree.RunNumber >= 188902]
            mc_eff = abs(corr.get3DMCEff(
                tau.pt, tau.eta,
                npileup_vtx, 0))
            mc_eff_errup = abs(corr.get3DMCEff(
                tau.pt, tau.eta,
                npileup_vtx, 1))
            mc_eff_errdn = abs(corr.get3DMCEff(
                tau.pt, tau.eta,
                npileup_vtx, -1))

            sf = corr.getSF(tau.pt, 0)
            sf_errup = corr.getSF(tau.pt, 1)
            sf_errdn = corr.getSF(tau.pt, -1)

            eff = mc_eff * sf
            try:
                eff_high = eff + eff * sqrt((mc_eff_errup / mc_eff)**2 + (sf_errup / sf)**2)
                eff_low = eff - eff * sqrt((mc_eff_errdn / mc_eff)**2 + (sf_errdn / sf)**2)
            except ZeroDivisionError:
                eff_high = 0.
                eff_low = 0.

            tau.trigger_eff = eff
            tau.trigger_eff_high = eff_high
            tau.trigger_eff_low = eff_low
            tau.trigger_eff_stat_high = eff
            tau.trigger_eff_stat_low = eff
            tau.trigger_eff_sys_high = eff
            tau.trigger_eff_sys_low = eff

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
        assert len(event.taus) == 2
        assert event.taus[0].pt >= event.taus[1].pt

        period = self.get_period_2012(self.tree.RunNumber)
        thresh = []

        # new tool only accepts EVnone
        eveto = 'EVnone'

        # taus are already sorted in descending order by pT by TauLeadSublead
        for tau, ttc in zip(event.taus, (self.ttc_29, self.ttc_20)):

            # only correct matched taus
            if not tau.matched:
                continue

            wpflag = self.get_id_12(tau)

            if tau.numTrack > 1:
                prong = '3p'
            else:
                prong = '1p'

            # data efficiency
            eff = abs(ttc.getDataEff(tau.pt, tau.eta, 0, self.tree.RunNumber, prong, wpflag, eveto))

            if math.isinf(eff) or math.isnan(eff):
                log.warning("trigger data efficiency is infinite or NaN! Using 0.")
                tau.trigger_eff = 0.
                tau.trigger_eff_high = 0.
                tau.trigger_eff_low = 0.
                tau.trigger_eff_stat_high = 0.
                tau.trigger_eff_stat_low = 0.
                tau.trigger_eff_sys_high = 0.
                tau.trigger_eff_sys_low = 0.

            else:
                tau.trigger_eff = eff
                tau.trigger_eff_high = eff
                tau.trigger_eff_low = eff

                # Data stat uncert
                stat_up = abs(ttc.getDataEff(tau.pt, tau.eta, 1, self.tree.RunNumber, prong, wpflag, eveto))
                stat_dn = abs(ttc.getDataEff(tau.pt, tau.eta, -1, self.tree.RunNumber, prong, wpflag, eveto))
                tau.trigger_eff_stat_high = eff + stat_up
                tau.trigger_eff_stat_low = eff - stat_dn

                # Systematic uncert
                sys_up = abs(ttc.getDataEff(tau.pt, tau.eta, 2, self.tree.RunNumber, prong, wpflag, eveto))
                sys_dn = abs(ttc.getDataEff(tau.pt, tau.eta, -2, self.tree.RunNumber, prong, wpflag, eveto))
                tau.trigger_eff_sys_high = eff + sys_up
                tau.trigger_eff_sys_low = eff - sys_dn

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
                        "division by zero in trigger scale factors: using 0.")
                    sf = 0.
                if math.isinf(sf) or math.isnan(sf):
                    log.warning("trigger data efficiency is infinite or NaN! Using 0.")
                    tau.trigger_sf = 0.
                    tau.trigger_sf_high = 0.
                    tau.trigger_sf_low = 0.
                    tau.trigger_sf_mc_stat_high = 0.
                    tau.trigger_sf_mc_stat_low = 0.
                    tau.trigger_sf_data_stat_high = 0.
                    tau.trigger_sf_data_stat_low = 0.
                    tau.trigger_sf_sys_high = 0.
                    tau.trigger_sf_sys_low = 0.

                else:
                    tau.trigger_sf = sf
                    tau.trigger_sf_high = sf
                    tau.trigger_sf_low = sf
                    tau.trigger_sf_mc_stat_high = sf
                    tau.trigger_sf_mc_stat_low = sf
                    tau.trigger_sf_data_stat_high = sf
                    tau.trigger_sf_data_stat_low = sf
                    tau.trigger_sf_sys_high = sf
                    tau.trigger_sf_sys_low = sf

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
                    tau.trigger_sf = 0.
                    tau.trigger_sf_high = 0.
                    tau.trigger_sf_low = 0.
                    tau.trigger_sf_mc_stat_high = 0.
                    tau.trigger_sf_mc_stat_low = 0.
                    tau.trigger_sf_data_stat_high = 0.
                    tau.trigger_sf_data_stat_low = 0.
                    tau.trigger_sf_sys_high = 0.
                    tau.trigger_sf_sys_low = 0.

                else:
                    tau.trigger_sf = sf
                    tau.trigger_sf_high = sf
                    tau.trigger_sf_low = sf

                    # Data stat uncert
                    data_stat_up = abs(ttc.getSF(tau.pt, tau.eta, 1, period, prong, wpflag, eveto))
                    data_stat_dn = abs(ttc.getSF(tau.pt, tau.eta, -1, period, prong, wpflag, eveto))
                    tau.trigger_sf_data_stat_high = sf + data_stat_up
                    tau.trigger_sf_data_stat_low = sf - data_stat_dn

                    # MC stat uncert
                    mc_stat_up = abs(ttc.getSF(tau.pt, tau.eta, 2, period, prong, wpflag, eveto))
                    mc_stat_dn = abs(ttc.getSF(tau.pt, tau.eta, -2, period, prong, wpflag, eveto))
                    tau.trigger_sf_mc_stat_high = sf + mc_stat_up
                    tau.trigger_sf_mc_stat_low = sf - mc_stat_dn

                    # Systematic uncert
                    sys_up = abs(ttc.getSF(tau.pt, tau.eta, 3, period, prong, wpflag, eveto))
                    sys_dn = abs(ttc.getSF(tau.pt, tau.eta, -3, period, prong, wpflag, eveto))
                    tau.trigger_sf_sys_high = sf + sys_up
                    tau.trigger_sf_sys_low = sf - sys_dn

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

            eff = abs(ttc.getDataEff(tau.pt, tau.eta, 0, self.tree.RunNumber, prong, wpflag, eveto))

            if math.isinf(eff) or math.isnan(eff):
                log.warning("trigger data efficiency is infinite or NaN! Using 0.")
                tau.trigger_eff = 0.
                tau.trigger_sf_high = 0.
                tau.trigger_sf_low = 0.
                tau.trigger_sf_mc_stat_high = 0.
                tau.trigger_sf_mc_stat_low = 0.
                tau.trigger_sf_data_stat_high = 0.
                tau.trigger_sf_data_stat_low = 0.
                tau.trigger_sf_sys_high = 0.
                tau.trigger_sf_sys_low = 0.

            else:
                tau.trigger_eff = eff
                tau.trigger_eff_high = eff
                tau.trigger_eff_low = eff

                # Data stat uncert
                stat_up = abs(ttc.getDataEff(tau.pt, tau.eta, 1, self.tree.RunNumber, prong, wpflag, eveto))
                stat_dn = abs(ttc.getDataEff(tau.pt, tau.eta, -1, self.tree.RunNumber, prong, wpflag, eveto))
                tau.trigger_eff_stat_high = eff + stat_up
                tau.trigger_eff_stat_low = eff - stat_dn

                # Systematic uncert
                sys_up = abs(ttc.getDataEff(tau.pt, tau.eta, 2, self.tree.RunNumber, prong, wpflag, eveto))
                sys_dn = abs(ttc.getDataEff(tau.pt, tau.eta, -2, self.tree.RunNumber, prong, wpflag, eveto))
                tau.trigger_eff_sys_high = eff + sys_up
                tau.trigger_eff_sys_low = eff - sys_dn

        return True

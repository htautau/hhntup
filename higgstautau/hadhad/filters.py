from rootpy.tree.filtering import EventFilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
import math
from .models import TrueTauBlock
from . import track_counting
from .. import tauid
from ..tauid import IDLOOSE, IDMEDIUM, IDTIGHT
from . import log; log = log[__name__]


class TauIDSelection(EventFilter):

    def __init__(self, year, tree, **kwargs):

        self.tree = tree
        if year == 2011:
            self.passes = self.passes_2011
        elif year == 2012:
            self.passes = self.passes_2012
        else:
            raise ValueError('no tau ID selection defined for year %d' % year)
        super(TauIDSelection, self).__init__(**kwargs)

    def passes_2011(self, event):

        tau1, tau2 = event.taus
        if not (tau1.JetBDTSigLoose and tau2.JetBDTSigLoose):
            return False
        # signal region is: both taus medium
        # need loose taus for OSFF model
        if tau1.JetBDTSigMedium:
            tau1.id = IDMEDIUM
        else:
            tau1.id = IDLOOSE
        if tau2.JetBDTSigMedium:
            tau2.id = IDMEDIUM
        else:
            tau2.id = IDLOOSE
        self.tree.taus_pass = tau1.JetBDTSigMedium and tau2.JetBDTSigMedium
        return True

    def passes_2012(self, event):

        # taus must be sorted in descending order by pT
        tau1, tau2 = event.taus
        if not (tau1.JetBDTSigLoose and tau2.JetBDTSigLoose):
            return False
        # signal region is: both medium with at least one being tight
        if ((tau1.JetBDTSigMedium and tau2.JetBDTSigMedium) and
            (tau1.JetBDTSigTight or tau2.JetBDTSigTight)):
            # if both are tight then assign medium to one at random
            # so we can apply the SFs in an inclusive manner
            if tau1.JetBDTSigTight and tau2.JetBDTSigTight:
                if event.EventNumber % 2 == 1: # ODD
                    tau1.id = IDTIGHT
                    tau2.id = IDMEDIUM
                else: # EVEN
                    tau1.id = IDMEDIUM
                    tau2.id = IDTIGHT
            elif tau1.JetBDTSigTight:
                tau1.id = IDTIGHT
                tau2.id = IDMEDIUM
            else:
                tau1.id = IDMEDIUM
                tau2.id = IDTIGHT
            self.tree.taus_pass = True
        else:
            # need loose taus for OSFF model
            tau1.id = IDLOOSE
            tau2.id = IDLOOSE
            self.tree.taus_pass = False
        return True


class TauLeadSublead(EventFilter):

    def __init__(self, lead=35*GeV, sublead=25*GeV, **kwargs):

        super(TauLeadSublead, self).__init__(**kwargs)
        # Leading and subleading tau pT thresholds
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

        # TODO use tau27Ti_m1_tau18Ti_m1_L2loose for period E
        # need emulaion, SFs for this


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


class TaudR(EventFilter):

    def __init__(self, dr=3.2, **kwargs):

        super(TaudR, self).__init__(**kwargs)
        self.dr = dr

    def passes(self, event):

        assert len(event.taus) == 2
        tau1, tau2 = event.taus
        return utils.dR(tau1.eta, tau1.phi, tau2.eta, tau2.phi) < self.dr


class TruthMatching(EventFilter):

    def passes(self, event):

        for tau in event.taus:
            tau.matched = tau.trueTauAssoc_index > -1
        return True


class TauTrackRecounting(EventFilter):

    def __init__(self, year, datatype, **kwargs):

        self.year = year
        self.datatype = datatype
        super(TauTrackRecounting, self).__init__(**kwargs)

    def passes(self, event):

        for tau in event.taus:
            tau.numTrack_recounted = track_counting.count_tracks(
                    tau, event, self.year, self.datatype)
        return True


class EfficiencyScaleFactors(EventFilter):

    def __init__(self, year, passthrough=False, **kwargs):

        if not passthrough:
            if year == 2011:
                log.info("will apply 2011 ID SFs")
                self.passes = self.passes_2011
            elif year == 2012:
                log.info("will apply 2012 ID SFs")
                from externaltools.bundle_2012 import TauCorrUncert as TCU
                from ROOT import TauCorrUncert
                self.tool = TauCorrUncert.TauSF(TCU.RESOURCE_PATH)
                self.tool_ns = TauCorrUncert
                self.passes = self.passes_2012
            else:
                raise ValueError("No efficiency SFs for year %d" % year)
        super(EfficiencyScaleFactors, self).__init__(
            passthrough=passthrough, **kwargs)

    def get_id_2011(self, tau):

        if tau.id == IDLOOSE:
            return 'loose'
        elif tau.id == IDMEDIUM:
            return 'medium'
        elif tau.id == IDTIGHT:
            return 'tight'
        raise ValueError("tau is not loose, medium, or tight")

    def get_id_2012(self, tau):

        if tau.id == IDLOOSE:
            return self.tool_ns.BDTLOOSE
        elif tau.id == IDMEDIUM:
            return self.tool_ns.BDTMEDIUM
        elif tau.id == IDTIGHT:
            return self.tool_ns.BDTTIGHT
        raise ValueError("tau is not loose, medium, or tight")

    def passes_2011(self, event):

        for tau in event.taus:

            if not tau.matched:
                continue

            wp = self.get_id_2011(tau)

            # efficiency scale factor
            effic_sf, err = tauid.effic_sf_uncert_exc(wp, tau, 2011)
            tau.efficiency_scale_factor =  effic_sf
            tau.efficiency_scale_factor_high = effic_sf + err
            tau.efficiency_scale_factor_low = effic_sf - err

        return True

    def passes_2012(self, event):

        for tau in event.taus:

            if not tau.matched:
                continue

            wp = self.get_id_2012(tau)

            sf = self.tool.GetIDSF(wp, tau.eta, tau.numTrack, tau.pt)
            sf_stat = self.tool.GetIDStatUnc(wp, tau.eta, tau.numTrack, tau.pt)
            sf_sys = self.tool.GetIDSysUnc(wp, tau.eta, tau.numTrack, tau.pt)

            sf_uncert = sqrt(sf_stat**2 + sf_sys**2)

            tau.efficiency_scale_factor =  sf
            tau.efficiency_scale_factor_high = sf + sf_uncert
            tau.efficiency_scale_factor_low = sf - sf_uncert

        return True


class FakeRateScaleFactors(EventFilter):

    def __init__(self, year, datatype, tes_up_systematic=False, tes_down_systematic=False, passthrough=False, **kwargs):

        self.tes_up = tes_up_systematic
        self.tes_down = tes_down_systematic
        if tes_up_systematic is None: self.tes_up = False
        if tes_down_systematic is None: self.tes_down = False

        if not passthrough:
            self.year = year % 1000
            self.datatype = datatype
            if self.year == 11:
                from externaltools.bundle_2011 import TauFakeRates
                from ROOT import TauFakeRates as TFR
                fakerate_table = TauFakeRates.get_resource(
                        'FakeRateScaleFactor.txt')
                self.fakerate_tool = TFR.FakeRateScaler(fakerate_table)
                self.passes = self.passes_2011
                log.info("will apply 2011 fake rate scale factors")
            elif self.year == 12:
                from externaltools.bundle_2012 import TauFakeRates
                from ROOT import TauFakeRates as TFR
                self.fakerate_tool = TFR.FakeRateScaler()
                self.fakerate_tool.initialise(TauFakeRates.RESOURCE_PATH)
                self.fakerate_ns = TFR
                self.passes = self.passes_2012
                log.info("will apply 2012 fake rate scale factors")
            else:
                raise ValueError("No fakerates defined for year %d" % year)

        super(FakeRateScaleFactors, self).__init__(
                passthrough=passthrough, **kwargs)

    def get_id_2011(self, tau):

        # 2011 fake rates are inclusive
        if tau.id == IDLOOSE:
            return 'Loose'
        elif tau.id == IDMEDIUM:
            return 'Medium'
        elif tau.id == IDTIGHT:
            return 'Tight'
        raise ValueError("tau is not loose, medium, or tight")

    def get_id_2012(self, tau):

        # 2012 fake rates are exclusive
        if tau.JetBDTSigTight:
            return self.fakerate_ns.TIGHT
        elif tau.JetBDTSigMedium:
            return self.fakerate_ns.MEDIUM
        elif tau.JetBDTSigLoose:
            return self.fakerate_ns.LOOSE
        raise ValueError("tau is not loose, medium, or tight")

    def passes_2011(self, event):

        if event.RunNumber >= 188902:
            trig = "EF_tau%dT_medium1"
        else:
            trig = "EF_tau%d_medium1"

        for tau in event.taus:

            # fakerate only applies to taus that don't match truth
            if tau.matched:
                continue

            wpflag = self.get_id_2011(tau)

            sf = self.fakerate_tool.getScaleFactor(
                    tau.pt, wpflag,
                    trig % tau.trigger_match_thresh)
            tau.fakerate_scale_factor = sf
            tau.fakerate_scale_factor_high = (sf +
                    self.fakerate_tool.getScaleFactorUncertainty(
                        tau.pt, wpflag,
                        trig % tau.trigger_match_thresh, True))
            tau.fakerate_scale_factor_low = (sf -
                    self.fakerate_tool.getScaleFactorUncertainty(
                        tau.pt, wpflag,
                        trig % tau.trigger_match_thresh, False))
        return True

    def passes_2012(self, event):

        for tau in event.taus:

            # fakerate only applies to taus that don't match truth
            if tau.matched:
                continue

            # Get the reco SF
            sf_reco = self.fakerate_tool.getRecoSF(
                tau.pt, tau.numTrack, event.RunNumber)

            tau.fakerate_scale_factor_reco = sf_reco
            tau.fakerate_scale_factor_reco_high = sf_reco
            tau.fakerate_scale_factor_reco_low = sf_reco

            wpflag = self.get_id_2012(tau)

            tes_up = self.tes_up
            tes_down = self.tes_down

            if tau.trigger_match_thresh == 20:
                trigger = self.fakerate_ns.TAU20Ti
            elif tau.trigger_match_thresh == 29:
                trigger = self.fakerate_ns.TAU29Ti
            else:
                raise ValueError("trigger threshold %d not understood" %
                    tau.trigger_match_thresh)

            # using LOOSE lepton veto
            sf_numer = self.fakerate_tool.getEffData(
                tau.pt, tau.numTrack, event.RunNumber,
                wpflag, self.fakerate_ns.LOOSE, trigger)

            sf_numer_up = self.fakerate_tool.getEffDataUncertainty(
                tau.pt, tau.numTrack, event.RunNumber,
                wpflag, self.fakerate_ns.LOOSE, trigger, True)

            sf_numer_dn = self.fakerate_tool.getEffDataUncertainty(
                tau.pt, tau.numTrack, event.RunNumber,
                wpflag, self.fakerate_ns.LOOSE, trigger, False)

            if self.datatype == datasets.MC:
                sf_denom = self.fakerate_tool.getEffMC(
                    tau.pt, tau.numTrack, event.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE, trigger, tes_down, tes_up)

                sf_denom_up = self.fakerate_tool.getEffMCUncertainty(
                    tau.pt, tau.numTrack, event.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE, trigger, True, tes_down, tes_up)

                sf_denom_dn = self.fakerate_tool.getEffMCUncertainty(
                    tau.pt, tau.numTrack, event.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE, trigger, False, tes_down, tes_up)

            else: # embedding: no trigger in denominator
                sf_denom = self.fakerate_tool.getEffData(
                    tau.pt, tau.numTrack, event.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE,
                    self.fakerate_ns.TRIGGERNONE)

                sf_denom_up = self.fakerate_tool.getEffDataUncertainty(
                    tau.pt, tau.numTrack, event.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE,
                    self.fakerate_ns.TRIGGERNONE, True)

                sf_denom_dn = self.fakerate_tool.getEffDataUncertainty(
                    tau.pt, tau.numTrack, event.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE,
                    self.fakerate_ns.TRIGGERNONE, False)

            #log.info("data eff: %f, mc eff %f, wp %s, pt %f, ntrack %d, run: %d, trigger: %d" % (
            #    sf_data, sf_mc, wp, tau.pt, tau.numTrack, event.RunNumber,
            #    tau.trigger_match_thresh))

            if sf_numer == 0 or sf_denom == 0:
                log.warning("fake rate bug: efficiency == 0")
                sf = 1.
                sf_high = 1.
                sf_low = 1.
            else:
                sf = sf_numer / sf_denom

                sf_up = sf * math.sqrt(
                    (sf_denom_up / sf_denom)**2 +
                    (sf_numer_up / sf_numer)**2)

                sf_dn = sf * math.sqrt(
                    (sf_denom_dn / sf_denom)**2 +
                    (sf_numer_dn / sf_numer)**2)

                sf_high = sf + sf_up
                sf_low = sf - sf_dn

            if sf_low < 0:
                sf_low = 0.

            tau.fakerate_scale_factor = sf
            # uncertainty
            tau.fakerate_scale_factor_high = sf_high
            tau.fakerate_scale_factor_low = sf_low
            #log.info("sf: %f, high: %f, low: %f" % (sf, sf_high, sf_low))

        return True

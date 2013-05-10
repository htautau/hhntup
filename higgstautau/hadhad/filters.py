from rootpy.tree.filtering import EventFilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *
import math
from .models import TrueTauBlock
from . import track_counting
from .. import tauid
from . import log; log = log[__name__]


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


class Triggers(EventFilter):
    """
    See lowest unprescaled triggers here:
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled#Taus_electron_muon_MET
    """
    def __init__(self, year, old_skim=False, **kwargs):

        if year == 2011:
            if old_skim:
                self.passes = self.passes_11_old
            else:
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

    def passes_11_old(self, event):
        try:
            if 177986 <= event.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1_EMULATED
            elif 188902 <= event.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1_EMULATED
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

    def __init__(self, year, **kwargs):

        self.year = year
        super(EfficiencyScaleFactors, self).__init__(**kwargs)

    def passes(self, event):

        for tau in event.taus:
            if not tau.matched:
                continue

            wplevels = []
            wplevels.append('loose')
            if tau.JetBDTSigMedium:
                wplevels.append('medium')
            if tau.JetBDTSigTight:
                wplevels.append('tight')

            for wp in wplevels:
                # efficiency scale factor
                effic_sf, err = tauid.effic_sf_uncert_exc(wp, tau, self.year)
                setattr(tau, 'id_eff_sf_%s' % wp, effic_sf)
                # ALREADY ACCOUNTED FOR IN TauBDT SYSTEMATIC
                setattr(tau, 'id_eff_sf_%s_high' % wp, effic_sf + err)
                setattr(tau, 'id_eff_sf_%s_low' % wp, effic_sf - err)
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

    def passes_2011(self, event):

        if event.RunNumber >= 188902:
            trig = "EF_tau%dT_medium1"
        else:
            trig = "EF_tau%d_medium1"

        for tau in event.taus:

            # fakerate only applies to taus that don't match truth
            if tau.matched:
                continue

            wplevels = []
            wplevels.append('loose')
            if tau.JetBDTSigMedium:
                wplevels.append('medium')
            if tau.JetBDTSigTight:
                wplevels.append('tight')

            for wp in wplevels:
                wpflag = wp.capitalize()
                sf = self.fakerate_tool.getScaleFactor(
                        tau.pt, wpflag,
                        trig % tau.trigger_match_thresh)
                setattr(tau, 'fakerate_sf_%s' % wp, sf)
                setattr(tau, 'fakerate_sf_%s_high' % wp, (sf +
                        self.fakerate_tool.getScaleFactorUncertainty(
                            tau.pt, wpflag,
                            trig % tau.trigger_match_thresh, True)))
                setattr(tau, 'fakerate_sf_%s_low' % wp, (sf -
                        self.fakerate_tool.getScaleFactorUncertainty(
                            tau.pt, wpflag,
                            trig % tau.trigger_match_thresh, False)))
        return True

    def passes_2012(self, event):

        for tau in event.taus:

            # fakerate only applies to taus that don't match truth
            if tau.matched:
                continue

            # Get the reco SF
            sf_reco = self.fakerate_tool.getRecoSF(
                tau.pt, tau.numTrack, event.RunNumber)

            setattr(tau, 'fakerate_sf_reco', sf_reco)
            setattr(tau, 'fakerate_sf_reco_high', sf_reco)
            setattr(tau, 'fakerate_sf_reco_low', sf_reco)

            tes_up = self.tes_up
            tes_down = self.tes_down
            wplevels = dict()
            wplevels['loose'] = self.fakerate_ns.LOOSE
            if tau.JetBDTSigMedium:
                wplevels['medium'] = self.fakerate_ns.MEDIUM
            if tau.JetBDTSigTight:
                wplevels['tight'] = self.fakerate_ns.TIGHT

            if tau.trigger_match_thresh == 20:
                trigger = self.fakerate_ns.TAU20Ti
            elif tau.trigger_match_thresh == 29:
                trigger = self.fakerate_ns.TAU29Ti
            else:
                raise ValueError("trigger threshold %d not understood" %
                    tau.trigger_match_thresh)

            for wp, wpflag in wplevels.items():
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
                setattr(tau, 'fakerate_sf_%s' % wp, sf)

                # uncertainty
                setattr(tau, 'fakerate_sf_%s_high' % wp, sf_high)
                setattr(tau, 'fakerate_sf_%s_low' % wp, sf_low)

                #log.info("sf: %f, high: %f, low: %f" % (sf, sf_high, sf_low))
        return True

from rootpy.tree.filtering import EventFilter

from math import *
import math

from .. import utils
from ..units import GeV
from .. import datasets
from . import track_counting
from .. import tauid
from ..tauid import IDLOOSE, IDMEDIUM, IDTIGHT
from . import log; log = log[__name__]

class LeptonSelection(EventFilter):
    def __init__(self,min_leps,**kwargs):
        super(LeptonSelection,self).__init__(**kwargs)
        self.min_leps = min_leps
        #require lepton candidates
    def passes(self, event):
        return len(event.electrons) + len(event.muons) >=self.min_leps

class LepPT(EventFilter):
    def __init__(self, min_leps, el_thresh=15*GeV, mu_thresh=10 * GeV, **kwargs):
        self.min_leps = min_leps
        self.el_thresh = el_thresh
        self.mu_thresh = mu_thresh
        super(LepPT, self).__init__(**kwargs)
    def passes(self, event):
        event.electrons.select(lambda electron: electron.cl_E/cosh(electron.tracketa) > self.el_thresh)
        event.muons.select(lambda muon: muon.pt > self.mu_thresh)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class LepEta(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(LepEta, self).__init__(**kwargs)
    def passes(self, event):                                                                                         
        # both calo and leading track eta within 2.47                                                                
        event.electrons.select(lambda electron: 
            abs(electron.cl_eta) < 1.37 or
            1.52 < abs(electron.cl_eta) < 2.47)
        event.muons.select(lambda muon: abs(muon.eta) < 2.5)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class LepHad(EventFilter):
    def __init__(self, **kwargs):
        super(LepHad, self).__init__(**kwargs)
    def passes(self, event):
        #make sure it's only one had and one lep
        #sort in decending order by pT
        event.electrons.sort(key=lambda electron: electron.cl_E/cosh(electron.tracketa), reverse=True)
        event.muons.sort(key=lambda muon: muon.pt, reverse=True)
        event.taus.sort(key=lambda tau: tau.pt, reverse=True)
        if len(event.taus) != 1: 
            return False
        if len(event.electrons) == 1:
            return True
        elif len(event.muons) == 1:
            return True
        else: 
            return False

#from here on we get specific muon and electron filters.

class MuonID(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(MuonID, self).__init__(**kwargs)
    def passes(self, event):
        #misc muon cuts
        event.muons.select(lambda muon: muon.tight)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class MuonHits(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(MuonHits, self).__init__(**kwargs)
    def passes(self, event):
        # Bottom layer requirements 
        event.muons.select(lambda muon:
            muon.expectBLayerHit == 0 or muon.nBLHits > 0)
        # Pixel hits
        event.muons.select(lambda muon:
            (muon.nPixHits + muon.nPixelDeadSensors) > 0)
        # SCT hits
        event.muons.select(lambda muon:
            (muon.nSCTHits + muon.nSCTDeadSensors) > 0)
        # Holes
        event.muons.select(lambda muon:
            (muon.nPixHoles + muon.nSCTHoles) < 3)
        # TRT
        event.muons.select(lambda muon:
            (0.1 < abs(muon.eta) < 1.9 and
             (muon.nTRTHits+muon.nTRTOutliers) > 5 and
             muon.nTRTOutliers < 0.9 * (muon.nTRTHits+muon.nTRTOutliers)) or
             (( abs(muon.eta) <= 0.1 or abs(muon.eta) >= 1.9) and
             not ((muon.nTRTHits+muon.nTRTOutliers) > 5 and 
             muon.nTRTOutliers >= 0.9 * (muon.nTRTHits+muon.nTRTOutliers))))
        return len(event.electrons) + len(event.muons) >= self.min_leps

class MuonZ0pv(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(MuonZ0pv, self).__init__(**kwargs)
    def passes(self, event):
        event.muons.select(lambda muon: muon.trackz0pv < 10)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class MuonCone(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(MuonCone, self).__init__(**kwargs)
    def passes(self, event):
        event.muons.select(lambda muon: muon.ptcone40 / muon.pt < 0.09)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class ElectronAuthor(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(ElectronAuthor, self).__init__(**kwargs)
    def passes(self, event):
        #extra misc electron cuts
        event.electrons.select(lambda electron: electron.author in (1, 3))
        return len(event.electrons) + len(event.muons) >= self.min_leps

class ElectronPP(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(ElectronPP, self).__init__(**kwargs)
    def passes(self, event):
        event.electrons.select(lambda electron: electron.mediumPP)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class ElectronOQ(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(ElectronOQ, self).__init__(**kwargs)
    def passes(self, event):
        #extra misc electron cuts
        event.electrons.select(lambda electron: (electron.OQ & 1446)==0)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class ElectronCone(EventFilter):
    def __init__(self, min_leps, **kwargs):
        self.min_leps = min_leps
        super(ElectronCone, self).__init__(**kwargs)
    def passes(self, event):
        #percentage energy in cone. electron.cl_E/cosh(electron.tracketa) is electron pt
        event.electrons.select(lambda electron: electron.ptcone40 / (electron.cl_E / cosh(electron.tracketa)) <0.17)
        return len(event.electrons) + len(event.muons) >= self.min_leps

class Triggers(EventFilter):
    """
    See lowest unprescaled triggers here:
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/LowestUnprescaled#Taus_electron_muon_MET
    """
    def __init__(self, year, tree, datatype, passthrough=False, **kwargs):
        if (not passthrough) and datatype == datasets.EMBED:
            raise ValueError("Cannot apply trigger on embedding samples")
        if year == 2011:
            if datatype == datasets.MC:
                self.passes = self.passes_11_mc
            else:
                self.passes = self.passes_11_data
        elif year == 2012:
            if datatype == datasets.MC:
                self.passes = self.passes_12_mc
            else:
                self.passes = self.passes_12_data
        else:
            raise ValueError("No triggers defined for year %d" % year)
        self.tree = tree
        super(Triggers, self).__init__(passthrough=passthrough, **kwargs)

    def passes_11_mc(self, event):
        try:
            if 177986 <= self.tree.RunNumber <= 187815: # Periods B-K
                self.tree.trigger = event.EF_tau29_medium1_tau20_medium1 == 1
                return True
            elif 188902 <= self.tree.RunNumber <= 191933: # Periods L-M
                self.tree.trigger = event.EF_tau29T_medium1_tau20T_medium1 == 1
                return True
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (self.tree.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" %
                         self.tree.RunNumber)

    def passes_11_data(self, event):
        try:
            if 177986 <= self.tree.RunNumber <= 187815: # Periods B-K
                return event.EF_tau29_medium1_tau20_medium1
            elif 188902 <= self.tree.RunNumber <= 191933: # Periods L-M
                return event.EF_tau29T_medium1_tau20T_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (self.tree.RunNumber, e)
            raise e
        raise ValueError("No trigger condition defined for run %s" %
                         self.tree.RunNumber)

    def passes_12_mc(self, event):
        try:
            self.tree.trigger = event.EF_tau29Ti_medium1_tau20Ti_medium1 == 1
            return True
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (self.tree.RunNumber, e)
            raise e
        # TODO use tau27Ti_m1_tau18Ti_m1_L2loose for period E
        # need emulaion, SFs for this

    def passes_12_data(self, event):
        try:
            return event.EF_tau29Ti_medium1_tau20Ti_medium1
        except AttributeError, e:
            print "Missing trigger for run %i: %s" % (self.tree.RunNumber, e)
            raise e
        # TODO use tau27Ti_m1_tau18Ti_m1_L2loose for period E
        # need emulaion, SFs for this


class TauFakeRateScaleFactors(EventFilter):

    def __init__(self, year, datatype, tree,
                 tes_up=False, tes_down=False,
                 passthrough=False, **kwargs):
        self.tes_up = tes_up
        self.tes_down = tes_down
        log.info("TauFakeRateScaleFactors: TES UP {0}".format(tes_up))
        log.info("TauFakeRateScaleFactors: TES DOWN {0}".format(tes_down))
        if not passthrough:
            self.year = year % 1000
            self.datatype = datatype
            self.tree = tree
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
        super(TauFakeRateScaleFactors, self).__init__(
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
        assert len(event.taus) == 2
        assert event.taus[0].pt >= event.taus[1].pt

        if self.tree.RunNumber >= 188902:
            trig = "EF_tau%dT_medium1"
        else:
            trig = "EF_tau%d_medium1"

        for tau, thresh in zip(event.taus, (29, 20)):

            # fakerate only applies to taus that don't match truth
            if tau.matched:
                continue

            wpflag = self.get_id_2011(tau)

            sf = self.fakerate_tool.getScaleFactor(
                tau.pt, wpflag,
                trig % thresh)
            tau.fakerate_sf = sf
            tau.fakerate_sf_high = (sf +
                self.fakerate_tool.getScaleFactorUncertainty(
                    tau.pt, wpflag,
                    trig % thresh, True))
            tau.fakerate_sf_low = (sf -
                self.fakerate_tool.getScaleFactorUncertainty(
                    tau.pt, wpflag,
                    trig % thresh, False))
        return True

    def passes_2012(self, event):
        assert len(event.taus) == 2
        assert event.taus[0].pt >= event.taus[1].pt

        for tau, trigger in zip(event.taus, [self.fakerate_ns.TAU29Ti, self.fakerate_ns.TAU20Ti]):
            # fakerate only applies to taus that don't match truth
            if tau.matched:
                continue
            # Get the reco SF
            sf_reco = self.fakerate_tool.getRecoSF(
                tau.pt, tau.numTrack, self.tree.RunNumber)
            tau.fakerate_sf_reco = sf_reco
            # NOTE: no uncertainty on this SF?
            # KG: yes there are, but getRecoSF didn't get updated to access them.
            # It's small so we will ignore it
            tau.fakerate_sf_reco_high = sf_reco
            tau.fakerate_sf_reco_low = sf_reco

            wpflag = self.get_id_2012(tau)

            tes_up = self.tes_up
            tes_down = self.tes_down

            # using LOOSE lepton veto
            sf_numer = self.fakerate_tool.getEffData(
                tau.pt, tau.numTrack, self.tree.RunNumber,
                wpflag, self.fakerate_ns.LOOSE, trigger)

            sf_numer_up = self.fakerate_tool.getEffDataUncertainty(
                tau.pt, tau.numTrack, self.tree.RunNumber,
                wpflag, self.fakerate_ns.LOOSE, trigger, True)

            sf_numer_dn = self.fakerate_tool.getEffDataUncertainty(
                tau.pt, tau.numTrack, self.tree.RunNumber,
                wpflag, self.fakerate_ns.LOOSE, trigger, False)

            if self.datatype == datasets.MC:
                sf_denom = self.fakerate_tool.getEffMC(
                    tau.pt, tau.numTrack, self.tree.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE, trigger, tes_down, tes_up)

                sf_denom_up = self.fakerate_tool.getEffMCUncertainty(
                    tau.pt, tau.numTrack, self.tree.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE, trigger, True, tes_down, tes_up)

                sf_denom_dn = self.fakerate_tool.getEffMCUncertainty(
                    tau.pt, tau.numTrack, self.tree.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE, trigger, False, tes_down, tes_up)

            else: # embedding: no trigger in denominator
                sf_denom = self.fakerate_tool.getEffData(
                    tau.pt, tau.numTrack, self.tree.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE,
                    self.fakerate_ns.TRIGGERNONE)

                sf_denom_up = self.fakerate_tool.getEffDataUncertainty(
                    tau.pt, tau.numTrack, self.tree.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE,
                    self.fakerate_ns.TRIGGERNONE, True)

                sf_denom_dn = self.fakerate_tool.getEffDataUncertainty(
                    tau.pt, tau.numTrack, self.tree.RunNumber,
                    wpflag, self.fakerate_ns.LOOSE,
                    self.fakerate_ns.TRIGGERNONE, False)

            #log.info("data eff: %f, mc eff %f, wp %s, pt %f, ntrack %d, run: %d, trigger: %d" % (
            #    sf_data, sf_mc, wp, tau.pt, tau.numTrack, self.tree.RunNumber,
            #    tau.trigger_match_thresh))

            if sf_numer == 0 or sf_denom == 0:
                log.warning("fake rate bug: efficiency == 0, using sf of 0")
                sf = 0.
                sf_high = 0.
                sf_low = 0.

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

            tau.fakerate_sf = sf
            # uncertainty
            tau.fakerate_sf_high = sf_high
            tau.fakerate_sf_low = sf_low
            #log.info("sf: %f, high: %f, low: %f" % (sf, sf_high, sf_low))
        return True

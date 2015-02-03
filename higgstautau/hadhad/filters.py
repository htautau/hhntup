import ROOT
from rootpy.tree.filtering import EventFilter

from math import *
import math

from .models import TrueTauBlock
from .. import utils
from ..units import GeV
from .. import datasets
from . import track_counting
from .. import tauid
from ..tauid import IDLOOSE, IDMEDIUM, IDTIGHT
from . import log; log = log[__name__]


class TauIDSelection(EventFilter):

    def __init__(self, **kwargs):
        super(TauIDSelection, self).__init__(**kwargs)

    def passes(self, event):
        # Enum definition
        # https://svnweb.cern.ch/trac/atlasoff/browser/Event/xAOD/xAODTau/trunk/xAODTau/TauDefs.h#L96
        JetBDTSigMedium = ROOT.xAOD.TauJetParameters.JetBDTSigMedium
        JetBDTSigTight = ROOT.xAOD.TauJetParameters.JetBDTSigTight

        tau1, tau2 = event.taus
        # signal region is: both medium with at least one being tight
        if ((tau1.obj.isTau(JetBDTSigMedium) and tau2.obj.isTau(JetBDTSigMedium)) and
            (tau1.obj.isTau(JetBDTSigTight) or tau2.obj.isTau(JetBDTSigTight))):
            # if both are tight then assign medium to one at random
            # so we can apply the SFs in an inclusive manner
            if tau1.obj.isTau(JetBDTSigTight) and tau2.obj.isTau(JetBDTSigTight):
                if event.EventInfo.eventNumber() % 2 == 1: # ODD
                    tau1.id = IDTIGHT
                    tau2.id = IDMEDIUM
                else: # EVEN
                    tau1.id = IDMEDIUM
                    tau2.id = IDTIGHT
            elif tau1.obj.isTau(JetBDTSigTight):
                tau1.id = IDTIGHT
                tau2.id = IDMEDIUM
            else:
                tau1.id = IDMEDIUM
                tau2.id = IDTIGHT
            return True
        return False


class TauLeadSublead(EventFilter):

    def __init__(self, lead=35*GeV, sublead=25*GeV, **kwargs):
        super(TauLeadSublead, self).__init__(**kwargs)
        # Leading and subleading tau pT thresholds
        self.lead = lead
        self.sublead = sublead

    def passes(self, event):
        # sort in descending order by pT
        event.taus.sort(key=lambda tau: tau.obj.pt(), reverse=True)
        # only keep leading two taus
        event.taus.slice(0, 2)
        # Event passes if the highest pT tau is above the leading
        # pT threshold and the next subleading tau pT is above the subleading pT theshold
        return event.taus[0].obj.pt() > self.lead and event.taus[1].obj.pt() > self.sublead


class Triggers(EventFilter):
    # NOT CONVERTED TO XAOD YET

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


class ElectronVeto(EventFilter):

    def __init__(self, el_sel='Medium', **kwargs):
        self.el_sel = el_sel
        super(ElectronVeto, self).__init__(**kwargs)

    def passes(self, event):
        for el in event.electrons:
            tracketa = el.trackParticle().eta()
            pt = el.caloCluster().e() / cosh(tracketa)
            if pt <= 15 * GeV:
                continue
            if not ((abs(tracketa) < 1.37) or (1.52 < abs(tracketa) < 2.47)):
                continue
            if el.author() not in (1, 3):
                continue
            if not abs(el.charge()) == 1:
                continue
            if el.passSelection(self.el_sel) != 1:
                continue
            from array import array as carray
            mask = carray('I', [ROOT.xAOD.EgammaParameters.BADCLUSELECTRON])
            if not el.isGoodOQ(mask[0]):
                continue
            return False
        return True


from ..filters import muon_has_good_track

class MuonVeto(EventFilter):

    def __init__(self, **kwargs):
        super(MuonVeto, self).__init__(**kwargs)
        self.muon_tool = ROOT.CP.MuonSelectionTool('MuonSelectionTool')
        # use loose quality and |eta|>2.5
        self.muon_tool.setProperty('int')('MuQuality', 2)
        self.muon_tool.initialize()
    def passes(self, event):
        for muon in event.muons:
            if muon.pt() <= 10 * GeV:
                continue
            if not self.muon_tool.accept(muon):
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
        return tau1.dr(tau2) < self.dr
        # return utils.dR(
        #     tau1.obj.eta(), tau1.obj.phi(), tau2.obj.eta(), tau2.obj.phi()) < self.dr


class TauTrackRecounting(EventFilter):
    # NOT CONVERTED TO XAOD YET

    def __init__(self, year, use_ntup_value=False, **kwargs):
        super(TauTrackRecounting, self).__init__(**kwargs)
        if year > 2011:
            if use_ntup_value:
                log.info("TauTrackRecounting using value from NTUP")
                self.count_func = track_counting.count_tracks_2012_p1443
            else:
                self.count_func = track_counting.count_tracks_2012
        else:
            self.count_func = track_counting.count_tracks_2011

    def passes(self, event):
        count_func = self.count_func
        for tau in event.taus:
            tau.numTrack_recounted = count_func(tau, event)
        return True


class TauIDScaleFactors(EventFilter):

    def __init__(self, year, passthrough=False, **kwargs):
        if not passthrough:
            if year == 2011:
                raise RuntimeWarning('TauIDScaleFactor 2011 is deprecated !')
            elif year == 2012:
                log.info("will apply 2012 ID SFs")
                from ROOT.TauAnalysisTools import TauEfficiencyCorrectionsTool
                # Medium id tool
                self.tool_medium = TauEfficiencyCorrectionsTool('tool_SF_med')
                self.tool_medium.setProperty('int')('IDLevel', 3)
                self.tool_medium.printConfig()
                self.tool_medium.initialize()

                # Tight id tool
                self.tool_tight = TauEfficiencyCorrectionsTool('tool_SF_tight')
                self.tool_tight.printConfig()
                self.tool_tight.initialize()

                self.passes = self.passes_2012
            else:
                raise ValueError("No efficiency SFs for year %d" % year)
        super(TauIDScaleFactors, self).__init__(
            passthrough=passthrough, **kwargs)

    def get_id_2011(self, tau):
        raise ValueError("deprecated !")

    def get_id_tool(self, tau):
        if tau.id == IDLOOSE:
            raise RuntimeError('Need to instantiate the ID loose tool !')
        elif tau.id == IDMEDIUM:
            return self.tool_medium
        elif tau.id == IDTIGHT:
            return self.tool_tight
        raise ValueError("tau is not loose, medium, or tight")

    def passes_2011(self, event):
        raise RuntimeError('deprecated !')
        return True

    def passes_2012(self, event):
        for tau in event.taus:
            if not tau.matched:
                continue
            tool = self.get_id_tool(tau)
            tool.applyEfficiencyScaleFactor(tau.obj)
            sf = tau.obj.auxdataConst('double')('TauScaleFactorJetID')
            systs = tool.recommendedSystematics()
            sf_syst = {}
            for sys in systs:
                # Why do we have to use a systematic set ??
                s_set = ROOT.CP.SystematicSet()
                s_set.insert(sys)
                tool.applySystematicVariation(s_set)
                tool.applyEfficiencyScaleFactor(tau.obj)
                sf_syst[sys.name()] = tau.obj.auxdataConst('double')('TauScaleFactorJetID')
            tau.id_sf =  sf
            # tau.id_sf_high = sf
            # tau.id_sf_low = sf
            tau.id_sf_stat_high = sf_syst['TAUS_EFF_JETID_STAT__1up']
            tau.id_sf_stat_low = sf_syst['TAUS_EFF_JETID_STAT__1down']
            tau.id_sf_sys_high = sf_syst['TAUS_EFF_JETID_SYST__1up']
            tau.id_sf_sys_low = sf_syst['TAUS_EFF_JETID_SYST__1down']
        return True


class TauFakeRateScaleFactors(EventFilter):
    # NOT CONVERTED TO XAOD YET

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
                # from externaltools.bundle_2011 import TauFakeRates
                from ROOT import TauFakeRates as TFR
                fakerate_table = TauFakeRates.get_resource(
                        'FakeRateScaleFactor.txt')
                self.fakerate_tool = TFR.FakeRateScaler(fakerate_table)
                self.passes = self.passes_2011
                log.info("will apply 2011 fake rate scale factors")
            elif self.year == 12:
                # from externaltools.bundle_2012 import TauFakeRates
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

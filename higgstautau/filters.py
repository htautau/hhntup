"""
Event filters common to both hadhad and lephad go here
"""

from rootpy.tree.filtering import *
from itertools import ifilter
from atlastools import utils
from atlastools.units import GeV
from atlastools import datasets
from math import *

from . import jetcleaning


def primary_vertex_selection(vxp):

    return vxp.type == 1 and vxp.nTracks >= 4


def pileup_vertex_selection(vxp):

    return vxp.type == 3 and vxp.nTracks >= 2


def vertex_selection(vxp):
    """ Does the full primary and pileup vertex selection """
    return primary_vertex_selection(vxp) or pileup_vertex_selection(vxp)


class PriVertex(EventFilter):

    def passes(self, event):

        event.vertices.select(vertex_selection)
        return any(ifilter(primary_vertex_selection, event.vertices))


class CoreFlags(EventFilter):

    def passes(self, event):

        return (event.coreFlags & 0x40000) == 0


class JetCleaning(EventFilter):

    def __init__(self,
                 datatype,
                 year,
                 level=jetcleaning.LOOSER,
                 pt_thresh=20*GeV,
                 eta_max=4.5,
                 **kwargs):

        super(JetCleaning, self).__init__(**kwargs)
        self.year = year
        self.datatype = datatype
        self.level = level
        self.pt_thresh = pt_thresh
        self.eta_max = eta_max

    def passes(self, event):

        # using LC jets
        for jet in event.jets:

            if jet.pt <= self.pt_thresh or abs(jet.eta) >= self.eta_max: continue
            LArQmean = jet.AverageLArQF / 65535.0
            chf = jet.sumPtTrk / jet.pt
            if jetcleaning.is_bad(
                    level=self.level,
                    quality=jet.LArQuality,
                    NegE=jet.NegativeE,
                    emf=jet.emfrac,
                    hecf=jet.hecf,
                    time=jet.Timing,
                    fmax=jet.fracSamplingMax,
                    eta=jet.emscale_eta,
                    chf=chf,
                    HecQ=jet.HECQuality,
                    LArQmean=LArQmean):
                return False

        if self.datatype == datasets.DATA and self.year == 2012:
            # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HowToCleanJets2012
            # Hot Tile calorimeter in period B1 and B2
            if 202660 <= event.RunNumber <= 203027:
                # recommendation is to use EM jets
                for jet in event.jets_EM:
                    _etaphi28 = (
                        -0.2 < jet.eta < -0.1 and
                        2.65 < jet.phi < 2.75)
                    if jet.fracSamplingMax > 0.6 and jet.SamplingMax == 13 and _etaphi28:
                        return False
            # Bad FCAL response in periods C1-C8
            # not applied in the skim for now.
            # Need to also apply on MC and choose a random run number with the
            # PileupReweighting tool
            #if 206248 <= event.RunNumber <= 207332:
            #    for jet in event.jets_EM:
            #        if (jet.pt > 20 * GeV and
            #            abs(jet.eta) > 3.2 and
            #            1.6 < jet.phi < 3.1):
            #            return False
        return True


class LArError(EventFilter):

    def passes(self, event):

        return event.larError <= 1


def in_lar_hole(eta, phi):

    return (-0.2 < eta < 1.6) and (-0.988 < phi < -0.392)


class LArHole(EventFilter):

    def __init__(self, datatype, **kwargs):

        super(LArHole, self).__init__(**kwargs)
        if datatype in (datasets.DATA, datasets.EMBED):
            self.passes = self.passes_data
        else:
            self.passes = self.passes_mc

    def passes_data(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        for jet in event.jets:
            if not jet.pt > 20 * GeV * (1 - jet.BCH_CORR_JET) / (1 - jet.BCH_CORR_CELL): continue
            if in_lar_hole(jet.eta, jet.phi): return False
        return True

    def passes_mc(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        for jet in event.jets:
            if not jet.pt > 20 * GeV: continue
            if in_lar_hole(jet.eta, jet.phi): return False
        return True


class JetCrackVeto(EventFilter):

    def passes(self, event):

        for jet in event.jets:
            if jet.pt <= 20 * GeV: continue
            if 1.3 < abs(jet.eta) < 1.7: return False
        return True


def muon_has_good_track(muon, year=2011):

    if year == 2011:
        pix_min = 2
        sct_min = 6
        abs_eta_min = -0.1
    elif year == 2012:
        pix_min = 1
        sct_min = 5
        abs_eta_min = 0.1
    else:
        raise ValueError("No muon veto defined for year %d" % year)

    blayer = (muon.expectBLayerHit == 0) or (muon.nBLHits > 0)
    pix = muon.nPixHits + muon.nPixelDeadSensors >= pix_min
    sct = muon.nSCTHits + muon.nSCTDeadSensors >= sct_min
    holes = muon.nPixHoles + muon.nSCTHoles < 3
    n_trt_hits_outliers = muon.nTRTHits + muon.nTRTOutliers

    if abs_eta_min < abs(muon.eta) < 1.9:
        trt = ((n_trt_hits_outliers > 5) and
              (muon.nTRTOutliers < 0.9 * n_trt_hits_outliers))
    else:
        trt = (n_trt_hits_outliers <= 5 or
               muon.nTRTOutliers < 0.9 * n_trt_hits_outliers)

    return blayer and pix and sct and holes and trt


class TauElectronVeto(EventFilter):

    def __init__(self, min_taus, **kwargs):

        super(TauElectronVeto, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):

        event.taus.select(lambda tau: tau.EleBDTLoose == 0)
        return len(event.taus) >= self.min_taus


class TauMuonVeto(EventFilter):

    def __init__(self, min_taus, **kwargs):

        super(TauMuonVeto, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):

        event.taus.select(lambda tau: tau.muonVeto == 0)
        return len(event.taus) >= self.min_taus


class TauHasTrack(EventFilter):

    def __init__(self, min_taus, **kwargs):

        super(TauHasTrack, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):

        event.taus.select(lambda tau: tau.numTrack > 0)
        return len(event.taus) >= self.min_taus


class TauAuthor(EventFilter):

    def __init__(self, min_taus, **kwargs):

        super(TauAuthor, self).__init__(**kwargs)
        self.min_taus = min_taus

    def passes(self, event):

        event.taus.select(lambda tau: tau.author != 2)
        return len(event.taus) >= self.min_taus


class TauPT(EventFilter):

    def __init__(self, min_taus, thresh=20 * GeV, **kwargs):

        self.min_taus = min_taus
        self.thresh = thresh
        super(TauPT, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: tau.pt > self.thresh)
        return len(event.taus) >= self.min_taus


class TauEta(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauEta, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: abs(tau.eta) < 2.5) # was 2.1
        return len(event.taus) >= self.min_taus


class TauJVF(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauJVF, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: True if
                abs(tau.track_eta[tau.leadtrack_idx]) > 2.1 else tau.jet_jvtxf > .5)
        return len(event.taus) >= self.min_taus


class Tau1Track3Track(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(Tau1Track3Track, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: tau.numTrack in (1, 3))
        return len(event.taus) >= self.min_taus


class TauCharge(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauCharge, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: abs(tau.charge) == 1)
        return len(event.taus) >= self.min_taus


class TauIDLoose(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauIDLoose, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: tau.JetBDTSigLoose)
        return len(event.taus) >= self.min_taus


class TauID_BDTLoose_LLHLoose(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauID_BDTLoose_LLHLoose, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau:
                tau.tauLlhLoose == 1 or tau.JetBDTSigLoose == 1)
        return len(event.taus) >= self.min_taus


class TauIDMedium(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauIDMedium, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: tau.JetBDTSigMedium)
        return len(event.taus) >= self.min_taus


class TauCrack(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauCrack, self).__init__(**kwargs)

    def passes(self, event):

        event.taus.select(lambda tau: not (1.37 < abs(tau.track_eta[tau.leadtrack_idx]) < 1.52))
        return len(event.taus) >= self.min_taus


class TauLArHole(EventFilter):

    def __init__(self, min_taus, **kwargs):

        self.min_taus = min_taus
        super(TauLArHole, self).__init__(**kwargs)

    def passes(self, event):

        if not 180614 <= event.RunNumber <= 184169:
            return True

        event.taus.select(lambda tau:
                not (-0.1 < tau.track_eta[tau.leadtrack_idx] < 1.55
                 and -0.9 < tau.track_phi[tau.leadtrack_idx] < -0.5))
        return len(event.taus) >= self.min_taus


def jet_selection(jet):
    """ Finalizes the jet selection """

    if not (jet.pt > 25*GeV) : return False

    #Protection against bunny ear jets
    if (2.5 < abs(jet.eta) < 3.5):
        if not (jet.pt > 30*GeV) : return False

    if not (abs(jet.eta) < 4.5) : return False
    if (abs(jet.eta) < 2.4):
        if not (jet.jvtxf > 0.75) : return False

    return True


class JetSelection(EventFilter):
    """Selects jets of good quality, keep event in any case"""

    def passes(self, event):

        event.jets.select(jet_selection)
        return True

"""
Event filters common to both hadhad and lephad go here
"""
import ROOT

from rootpy.tree.filtering import *
from itertools import ifilter
from math import *

from . import datasets
from .corrections import reweight_ggf
from .units import GeV
from . import jetcleaning
from . import utils
from . import log; log = log[__name__]

from goodruns import GRL


class GRLFilter(EventFilter):

    def __init__(self, grl, **kwargs):
        super(GRLFilter, self).__init__(**kwargs)
        if isinstance(grl, GRL):
            self.grl = grl
        else:
            self.grl = GRL(grl)

    def passes(self, event):
        if not self.grl:
            return False
        return (event.RunNumber, event.lbn) in self.grl


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


class RandomSeed(EventFilter):

    def __init__(self, datatype, **kwargs):
        super(RandomSeed, self).__init__(**kwargs)
        self.datatype = datatype

    def passes(self, event):
        # ResoSoftTerms uses gRandom for smearing.
        # Set the seed here however you like.
        if self.datatype in (datasets.DATA, datasets.EMBED):
            ROOT.gRandom.SetSeed(int(event.RunNumber * event.EventNumber))
        else:
            ROOT.gRandom.SetSeed(int(event.mc_channel_number * event.EventNumber))
        return True


class RandomRunNumber(EventFilter):

    def __init__(self, tree, datatype, pileup_tool, **kwargs):
        self.tree = tree
        self.pileup_tool = pileup_tool
        super(RandomRunNumber, self).__init__(**kwargs)
        if datatype in (datasets.MC, datasets.MCEMBED):
            self.passes = self.passes_mc
        else:
            self.passes = self.passes_data

    def passes_data(self, event):
        self.tree.RunNumber = event.RunNumber
        return True

    def passes_mc(self, event):
        # get random run number using the pileup tool
        self.tree.RunNumber = self.pileup_tool.GetRandomRunNumber(event.RunNumber)
        return True


class NvtxJets(EventFilter):

    def __init__(self, tree, **kwargs):
        super(NvtxJets, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        # Check for a good primary vertex
        # This is needed for jet and soft term systematics
        goodPV = False
        nvtxsoftmet = 0
        nvtxjets = 0
        # Most D3PDs contain the vx_type branch, but some don't.
        # Those which don't are most likely skimmed, and require at least 1
        # primary vertex for all events.
        # If your D3PD is skimmed in this way, then the goodPV (nTracks and z)
        # check should be applied to the first vertex in the collection.
        # Otherwise, you should ensure that the vx_type branch is available.
        for vertex in event.vertices:
            if vertex.type == 1 and vertex.nTracks > 2 and abs(vertex.z) < 200:
                goodPV = True
        if goodPV:
            for vertex in event.vertices:
                if vertex.nTracks > 2:
                    nvtxsoftmet += 1
                if vertex.nTracks > 1:
                    nvtxjets += 1
        self.tree.nvtxsoftmet = nvtxsoftmet
        self.tree.nvtxjets = nvtxjets
        return True


class TileTrips(EventFilter):
    """
    https://twiki.cern.ch/twiki/bin/viewauth/Atlas/DataPreparationCheckListForPhysicsAnalysis#Rejection_of_bad_corrupted_event
    """
    def __init__(self, passthrough=False, **kwargs):
        if not passthrough:
            from externaltools import TileTripReader
            from ROOT import Root
            self.tool = Root.TTileTripReader()
        super(TileTrips, self).__init__(passthrough=passthrough, **kwargs)

    def passes(self, event):
        # only apply between G - J
        #if event.RunNumber < 211522:
        #    return True
        #if event.RunNumber > 215091:
        #    return True
        # returns false if the event is one with a saturation in a tile cell
        # (bad MET).
        return self.tool.checkEvent(
            event.RunNumber,
            event.lbn,
            event.EventNumber)


class JetCleaning(EventFilter):

    # TODO:
    # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/BCHCleaningTool

    BAD_TILE = [
        202660, 202668, 202712, 202740, 202965, 202987, 202991, 203027, 203169
    ]

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
            if jet.pt <= self.pt_thresh or abs(jet.eta) >= self.eta_max:
                continue
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

        if (self.datatype in (datasets.DATA, datasets.EMBED)) and self.year == 2012:
            # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HowToCleanJets2012
            # Hot Tile calorimeter in period B1 and B2
            if event.RunNumber in JetCleaning.BAD_TILE:
                # recommendation is to use EM jets
                for jet in event.jets_EM:
                    _etaphi28 = (
                        -0.2 < jet.eta < -0.1 and
                        2.65 < jet.phi < 2.75)
                    if jet.fracSamplingMax > 0.6 and jet.SamplingMax == 13 and _etaphi28:
                        return False
            # Not required in reprocessed data:
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
        return event.larError != 2


class TileError(EventFilter):

    def passes(self, event):
        return event.tileError != 2


def in_lar_hole(eta, phi):
    return (-0.2 < eta < 1.6) and (-0.988 < phi < -0.392)


class LArHole(EventFilter):
    """
    https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HowToCleanJets2011#LAr_Hole
    """
    def __init__(self, tree, **kwargs):
        super(LArHole, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        # only apply from period E to H
        if not 180614 <= self.tree.RunNumber <= 184169:
            return True
        for jet in event.jets:
            if not jet.pt * (1. - jet.BCH_CORR_CELL) / (1. - jet.BCH_CORR_JET) > 20 * GeV:
                continue
            if in_lar_hole(jet.eta, jet.phi):
                return False
        return True


class JetCrackVeto(EventFilter):

    def passes(self, event):
        for jet in event.jets:
            if jet.pt <= 20 * GeV: continue
            if 1.3 < abs(jet.eta) < 1.7: return False
        return True


def muon_has_good_track(muon, year):
    """
    https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MCPAnalysisGuidelinesRel17MC11a
    https://twiki.cern.ch/twiki/bin/view/AtlasProtected/MCPAnalysisGuidelinesData2012
    """
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
        #https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauRecommendationsWinterConf2013#Electron_veto
        # only apply eveto on 1p taus with cluster and track eta less than 2.47
        # Eta selection already applied by TauEta filter
        event.taus.select(lambda tau:
            tau.numTrack > 1 or
            tau.EleBDTLoose == 0)
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
        # both calo and leading track eta within 2.47
        event.taus.select(lambda tau:
            abs(tau.eta) < 2.47 and
            abs(tau.track_eta[tau.leadtrack_idx]) < 2.47)
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


class Tau1P3P(EventFilter):
    """
    Only keep 1P + 3P and 3P + 3P
    """
    def passes(self, event):
        assert len(event.taus) == 2
        tau1, tau2 = event.taus
        # 1P + 3P
        if (tau1.numTrack == 1 and tau2.numTrack == 3) or \
           (tau1.numTrack == 3 and tau2.numTrack == 1):
            return True
        # 3P + 3P
        if tau1.numTrack == 3 and tau2.numTrack == 3:
            return True
        return False


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
        event.taus.select(lambda tau: tau.JetBDTSigLoose == 1)
        return len(event.taus) >= self.min_taus


class TauIDMedium(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauIDMedium, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: tau.JetBDTSigMedium == 1)
        return len(event.taus) >= self.min_taus


class TauCrack(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauCrack, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(
            lambda tau: not (
                1.37 <= abs(tau.track_eta[tau.leadtrack_idx]) <= 1.52))
        return len(event.taus) >= self.min_taus


class TauLArHole(EventFilter):

    def __init__(self, min_taus, tree, **kwargs):
        self.min_taus = min_taus
        self.tree = tree
        super(TauLArHole, self).__init__(**kwargs)

    def passes(self, event):
        if not 180614 <= self.tree.RunNumber <= 184169:
            return True
        event.taus.select(lambda tau:
            not (-0.1 < tau.track_eta[tau.leadtrack_idx] < 1.55
                 and -0.9 < tau.track_phi[tau.leadtrack_idx] < -0.5))
        return len(event.taus) >= self.min_taus


class TruthMatching(EventFilter):

    def passes(self, event):
        for tau in event.taus:
            # CANNOT do the following due to buggy D3PD:
            # tau.matched = tau.trueTauAssoc_index > -1
            tau.matched = False
            tau.matched_dr = 1111.
            tau.matched_object = None
            for truetau in event.truetaus:
                dr = utils.dR(tau.eta, tau.phi, truetau.vis_eta, truetau.vis_phi)
                if dr < 0.2:
                    # TODO: handle possible collision!
                    tau.matched = True
                    tau.matched_dr = dr
                    tau.matched_object = truetau
                    break
        return True


class TauSelected(EventFilter):

    def __init__(self, min_taus, **kwargs):
        self.min_taus = min_taus
        super(TauSelected, self).__init__(**kwargs)

    def passes(self, event):
        event.taus.select(lambda tau: tau.selected)
        return len(event.taus) >= self.min_taus


class TauEnergyShift(EventFilter):
    """
    in situ TES shift for 8TeV 2012 data
    """
    def __init__(self, *args, **kwargs):
        from externaltools import TauCorrUncert as TCU
        from ROOT import TauCorrUncert
        self.tool = TauCorrUncert.TESUncertainty(
            TCU.get_resource('TES/mc12_p1344_medium.root'))
        super(TauEnergyShift, self).__init__(*args, **kwargs)

    def passes(self, event):
        shift_func = self.tool.GetTESShift
        for tau in event.taus:
            shift = shift_func(tau.pt, tau.numTrack)
            tau.pt *= 1. + shift
        return True


class NumJets25(EventFilter):

    def __init__(self, tree, **kwargs):
        super(NumJets25, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        self.tree.numJets25 = len([j for j in event.jets if
            j.pt > 25 * GeV and abs(j.eta) < 4.5])
        return True


class NonIsolatedJet(EventFilter):
    """
    https://indico.cern.ch/getFile.py/access?contribId=1&resId=0&materialId=slides&confId=200403
    """
    def __init__(self, tree, **kwargs):
        super(NonIsolatedJet, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        # only write flag instead of vetoing the event so this
        # can be turned on and off after
        self.tree.nonisolatedjet = False
        for tau in event.taus:
            for jet in event.jets:
                if 0.4 < utils.dR(tau.eta, tau.phi, jet.eta, jet.phi) < 1.0:
                    self.tree.nonisolatedjet = True
        return True


def jet_selection_2011(jet):
    """ Finalizes the jet selection """

    if not (jet.pt > 25 * GeV):
        return False

    if not (abs(jet.eta) < 4.5):
        return False

    # suppress forward jets
    if (abs(jet.eta) > 2.4) and not (jet.pt > 30 * GeV):
        return False

    # JVF cut on central jets
    #if (abs(jet.eta) < 2.4) and not (abs(jet.jvtxf) > 0.75):
    #    return False
    # NO JVFUncertaintyTool for 2011!

    return True


def jet_selection_2012(jet):
    """ Finalizes the jet selection
    https://cds.cern.ch/record/1472547/files/ATL-COM-PHYS-2012-1202.pdf
    """
    if not (jet.pt > 30 * GeV):
        return False

    if not (abs(jet.eta) < 4.5):
        return False

    # suppress forward jets
    if (abs(jet.eta) > 2.4) and not (jet.pt > 35 * GeV):
        return False

    # JVF cut on central jets below 50 GeV
    if (jet.pt < 50 * GeV) and (abs(jet.constscale_eta) < 2.4):
        if not (abs(jet.jvtxf) > 0.5):
            return False

    return True


class JetSelection(EventFilter):
    """Selects jets of good quality, keep event in any case"""

    def __init__(self, year, **kwargs):
        if year == 2011:
            self.filter_func = jet_selection_2011
        elif year == 2012:
            self.filter_func = jet_selection_2012
        else:
            raise ValueError("No jet selection defined for year %d" % year)
        super(JetSelection, self).__init__(**kwargs)

    def passes(self, event):
        event.jets.select(self.filter_func)
        return True


class JetPreselection(EventFilter):

    def passes(self, event):
        event.jets.select(lambda jet: jet.pt > 20 * GeV)
        return True


class MCWeight(EventFilter):

    def __init__(self, datatype, tree, **kwargs):
        self.datatype = datatype
        self.tree = tree
        super(MCWeight, self).__init__(**kwargs)

    def passes(self, event):
        # set the event weights
        if self.datatype == datasets.MC:
            self.tree.mc_weight = event.mc_event_weight
        elif self.datatype == datasets.EMBED:
            # https://twiki.cern.ch/twiki/bin/viewauth/Atlas/EmbeddingTools
            # correct truth filter efficiency
            self.tree.mc_weight = event.mcevt_weight[0][0]
            # In 2011 mc_event_weight == mcevt_weight[0][0]
            # for embedding. But not so in 2012...
        return True


class ggFReweighting(EventFilter):

    def __init__(self, dsname, tree, **kwargs):
        self.dsname = dsname
        self.tree = tree
        super(ggFReweighting, self).__init__(**kwargs)

    def passes(self, event):
        self.tree.ggf_weight = reweight_ggf(event, self.dsname)
        return True


class EmbeddingCorrections(EventFilter):
    # https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/HiggsToTauTauEmbeddedCorrections2013
    def __init__(self, tree, passthrough=False, **kwargs):
        super(EmbeddingCorrections, self).__init__(passthrough=passthrough, **kwargs)
        if not passthrough:
            self.tree = tree
            from externaltools import EmbeddedCorrections
            from externaltools import TrigMuonEfficiency
            from externaltools import ElectronEfficiencyCorrection
            from externaltools import HSG4LepLepTriggerSF
            from externaltools import MuonEfficiencyCorrections
            self.tool = ROOT.EmbeddedCorrections.Embedded(
                EmbeddedCorrections.get_resource('2DMaps.root'),
                EmbeddedCorrections.RESOURCE_PATH,
                ElectronEfficiencyCorrection.RESOURCE_PATH,
                HSG4LepLepTriggerSF.RESOURCE_PATH,
                MuonEfficiencyCorrections.RESOURCE_PATH)

    def passes(self, event):
        self.tool.SetupEmbeddedEvent(
            event.mc_pt,
            event.mc_eta,
            event.mc_phi,
            event.mc_m,
            event.mc_pdgId,
            # possibly a random run number (for MC embedding)
            self.tree.RunNumber)
        # Retrieve the unfolding weight
        self.tree.embedding_reco_unfold = self.tool.GetEmbeddingRecoUnfolding()
        # Access the trigger unfolding weight
        self.tree.embedding_trigger_weight = self.tool.GetEmbeddingTriggerWeight()
        # Get the original mass of the dimuon event
        self.tree.embedding_dimuon_mass = self.tool.GetOriginalZ().M()
        return True


class JetCopy(EventFilter):

    def __init__(self, tree, **kwargs):
        super(JetCopy, self).__init__(**kwargs)
        self.tree = tree

    def passes(self, event):
        tree = self.tree
        tree.jet_E_original.clear()
        tree.jet_m_original.clear()
        tree.jet_pt_original.clear()
        tree.jet_eta_original.clear()
        tree.jet_phi_original.clear()
        for jet in event.jets:
            tree.jet_E_original.push_back(jet.E)
            tree.jet_m_original.push_back(jet.m)
            tree.jet_pt_original.push_back(jet.pt)
            tree.jet_eta_original.push_back(jet.eta)
            tree.jet_phi_original.push_back(jet.phi)
        return True
